# =============================================================================
# TCR Diversity Analysis: NanoString-inspired TRBV Proxy vs True Clonotype Diversity
# =============================================================================
#
# DESCRIPTION:
#   This module implements a NanoString-inspired TCR beta variable (TRBV) gene
#   usage diversity proxy and compares it to true clonotype-level diversity
#   derived from matched single-cell TCR sequencing.
#
#   The TRBV-based metric measures V-gene usage breadth/evenness (analagous to
#   the NanoString nCounter TCR Diversity approach). It CANNOT resolve CDR3
#   clonotypes, antigen specificity, or clonal expansion within a TRBV family.
#   Do not interpret it as "true TCR diversity."
#
# KEY FUNCTIONS:
#   shannon_diversity(counts)
#   normalized_shannon(counts)
#   simpson_diversity(counts)
#   clonality_from_shannon(counts)
#   compute_trbv_diversity(metadata, group_vars, min_cells, downsample_n)
#   compute_clonotype_diversity(metadata, group_vars, min_cells, downsample_n)
#   compare_diversity_metrics(trbv_results, clonotype_results)
#   benchmark_failure_modes(metadata, group_vars)
#   run_tcr_diversity_analysis(metadata, group_vars, min_cells, downsample_n)
#
# DATA LOADING HELPERS:
#   load_vdj_data(vdj_dir, sample_names, sample_order)
#   load_aggr_barcodes(aggr_dir, sample_names, sample_order)
#   build_tcr_metadata(vdj_dir, sample_names, sample_order,
#                      aggr_dir, additional_metadata)
#
# PLOTTING FUNCTIONS:
#   plot_trbv_stacked_bar(metadata, group_var, top_n_genes)
#   plot_clonotype_rank_abundance(metadata, group_var)
#   plot_h_trbv_vs_h_clonotype(diversity_table)
#   plot_clonality_comparison(diversity_table)
#   plot_trbv_heatmap(metadata, group_var, top_n_genes)
#   plot_within_trbv_diversity(benchmark_results)
#
# AUTHORS: Generated for the Sidiropoulos Lab / J1568 Batch 1 analysis
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(data.table)
  library(vegan)       # for diversity metrics (optional, fallback implemented)
  library(RColorBrewer)
  library(patchwork)   # for combining ggplots
  library(pheatmap)
})


# =============================================================================
# SECTION 1: CORE DIVERSITY METRIC FUNCTIONS
# =============================================================================

#' Compute Shannon diversity from a named count vector
#'
#' @param counts  Numeric vector of counts (NAs and zeros removed internally).
#' @param base    Logarithm base (default = natural log, base = exp(1)).
#' @return  Shannon diversity H >= 0. Returns NA if fewer than 2 non-zero counts.
shannon_diversity <- function(counts, base = exp(1)) {
  counts <- counts[!is.na(counts) & counts > 0]
  if (length(counts) < 2) return(NA_real_)
  p <- counts / sum(counts)
  -sum(p * log(p, base = base))
}

#' Compute normalized Shannon diversity (Pielou's evenness)
#'
#' H_norm = H / log(S), where S = number of detected categories.
#' Ranges [0, 1]; 1 = perfectly even, 0 = single dominant category.
#'
#' @param counts  Numeric vector of counts.
#' @param base    Logarithm base (default = natural log).
#' @return  Normalized Shannon (evenness) in [0, 1], or NA if S < 2.
normalized_shannon <- function(counts, base = exp(1)) {
  counts <- counts[!is.na(counts) & counts > 0]
  if (length(counts) < 2) return(NA_real_)
  H <- shannon_diversity(counts, base = base)
  H / log(length(counts), base = base)
}

#' Compute Simpson diversity (1 - D, probability that two random cells differ)
#'
#' @param counts  Numeric vector of counts.
#' @return  Simpson diversity in [0, 1]. Returns NA if fewer than 2 non-zero counts.
simpson_diversity <- function(counts) {
  counts <- counts[!is.na(counts) & counts > 0]
  if (length(counts) < 2) return(NA_real_)
  n <- sum(counts)
  if (n < 2) return(NA_real_)
  # Unbiased estimator: 1 - sum(n_i*(n_i-1)) / (n*(n-1))
  1 - sum(counts * (counts - 1)) / (n * (n - 1))
}

#' Compute clonality as 1 - normalized Shannon
#'
#' Clonality = 0 means perfectly even distribution; 1 means single dominant clone.
#' This is analogous to the metric used in Adaptive Immunosequencing outputs.
#'
#' @param counts  Numeric vector of counts.
#' @param base    Logarithm base (default = natural log).
#' @return  Clonality in [0, 1], or NA if S < 2.
clonality_from_shannon <- function(counts, base = exp(1)) {
  ev <- normalized_shannon(counts, base = base)
  if (is.na(ev)) return(NA_real_)
  1 - ev
}

#' Compute Gini coefficient from a count vector
#'
#' Gini = 0: perfect equality; Gini = 1: maximum inequality.
#'
#' @param counts  Numeric vector of counts.
#' @return  Gini coefficient in [0, 1].
gini_coefficient <- function(counts) {
  counts <- counts[!is.na(counts) & counts > 0]
  n <- length(counts)
  if (n < 2) return(NA_real_)
  counts <- sort(counts)
  2 * sum((seq_len(n)) * counts) / (n * sum(counts)) - (n + 1) / n
}

#' Downsample a metadata data frame to n_target cells per group
#'
#' @param metadata   Data frame with at minimum a column named by group_vars.
#' @param group_vars Character vector of grouping column names.
#' @param n_target   Target number of cells per group. Groups with fewer cells
#'                   are kept as-is.
#' @param seed       Random seed for reproducibility.
#' @return  Downsampled data frame.
downsample_cells <- function(metadata, group_vars, n_target, seed = 42) {
  set.seed(seed)
  # Use group_modify so nrow(.x) is available as a plain integer (not a
  # data-masking expression), avoiding the dplyr n() restriction in slice_sample().
  metadata %>%
    group_by(across(all_of(group_vars))) %>%
    group_modify(~ slice_sample(.x, n = min(n_target, nrow(.x)),
                                replace = FALSE)) %>%
    ungroup()
}


# =============================================================================
# SECTION 2: PER-GROUP DIVERSITY COMPUTATION
# =============================================================================

#' Compute TRBV usage diversity (NanoString-style proxy) per group
#'
#' IMPORTANT: This metric measures V-gene usage breadth/evenness.
#' It CANNOT resolve CDR3 clonotypes, antigen specificity, or clonal
#' expansion within a TRBV family. Do NOT interpret as "true TCR diversity."
#'
#' @param metadata    Data frame with columns: sample_id, trbv_gene, and any
#'                    columns listed in group_vars.
#' @param group_vars  Character vector of column names to group by (e.g.,
#'                    c("sample_id") or c("sample_id", "subset")).
#' @param min_cells   Minimum cells required per group (default 30).
#' @param downsample_n  If > 0, downsample each group to this many cells before
#'                      computing diversity (sensitivity analysis). Default 0 (off).
#' @return  Data frame with one row per group and TRBV diversity columns.
compute_trbv_diversity <- function(metadata,
                                   group_vars = "sample_id",
                                   min_cells  = 30,
                                   downsample_n = 0) {

  # Require trbv_gene column
  if (!"trbv_gene" %in% colnames(metadata)) {
    stop("metadata must contain column 'trbv_gene'")
  }

  # Optionally downsample
  if (downsample_n > 0) {
    metadata <- downsample_cells(metadata, group_vars, downsample_n)
  }

  metadata %>%
    # Filter rows with valid TRBV assignment
    filter(!is.na(trbv_gene), trbv_gene != "", trbv_gene != "None") %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      n_cells_trbv       = n(),
      n_trbv_detected    = n_distinct(trbv_gene),
      # TRBV gene counts → diversity metrics
      H_TRBV             = {
        cts <- table(trbv_gene)
        shannon_diversity(as.numeric(cts))
      },
      H_norm_TRBV        = {
        cts <- table(trbv_gene)
        normalized_shannon(as.numeric(cts))
      },
      TRBV_clonality_proxy = {
        cts <- table(trbv_gene)
        clonality_from_shannon(as.numeric(cts))
      },
      simpson_TRBV       = {
        cts <- table(trbv_gene)
        simpson_diversity(as.numeric(cts))
      },
      gini_TRBV          = {
        cts <- table(trbv_gene)
        gini_coefficient(as.numeric(cts))
      },
      top_TRBV_gene      = names(sort(table(trbv_gene), decreasing = TRUE))[1],
      top_TRBV_fraction  = {
        cts <- table(trbv_gene)
        max(cts) / sum(cts)
      },
      .groups = "drop"
    ) %>%
    # Apply minimum cell threshold
    filter(n_cells_trbv >= min_cells)
}

#' Compute true clonotype-level diversity per group
#'
#' Uses clonotype_id (full CDR3-based clonotype assignment) to compute
#' true clonal diversity. This is the "ground truth" metric.
#'
#' @param metadata    Data frame with columns: sample_id, clonotype_id, and any
#'                    columns listed in group_vars.
#' @param group_vars  Character vector of column names to group by.
#' @param min_cells   Minimum cells required per group (default 30).
#' @param downsample_n  If > 0, downsample each group before computing diversity.
#' @return  Data frame with one row per group and clonotype diversity columns.
compute_clonotype_diversity <- function(metadata,
                                        group_vars   = "sample_id",
                                        min_cells    = 30,
                                        downsample_n = 0) {

  if (!"clonotype_id" %in% colnames(metadata)) {
    stop("metadata must contain column 'clonotype_id'")
  }

  if (downsample_n > 0) {
    metadata <- downsample_cells(metadata, group_vars, downsample_n)
  }

  metadata %>%
    filter(!is.na(clonotype_id), clonotype_id != "", clonotype_id != "None") %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      n_cells_clonotype   = n(),
      n_clonotypes        = n_distinct(clonotype_id),
      H_clonotype         = {
        cts <- table(clonotype_id)
        shannon_diversity(as.numeric(cts))
      },
      H_norm_clonotype    = {
        cts <- table(clonotype_id)
        normalized_shannon(as.numeric(cts))
      },
      true_clonality      = {
        cts <- table(clonotype_id)
        clonality_from_shannon(as.numeric(cts))
      },
      simpson_clonotype   = {
        cts <- table(clonotype_id)
        simpson_diversity(as.numeric(cts))
      },
      gini_clonotype      = {
        cts <- table(clonotype_id)
        gini_coefficient(as.numeric(cts))
      },
      top_clone_id        = names(sort(table(clonotype_id), decreasing = TRUE))[1],
      top_clone_fraction  = {
        cts <- table(clonotype_id)
        max(cts) / sum(cts)
      },
      .groups = "drop"
    ) %>%
    filter(n_cells_clonotype >= min_cells)
}

#' Join TRBV proxy and clonotype diversity results and compute correlation
#'
#' @param trbv_results       Output of compute_trbv_diversity().
#' @param clonotype_results  Output of compute_clonotype_diversity().
#' @param group_vars         Character vector of join keys (must match both inputs).
#' @return  List with:
#'   $table       — joined diversity table
#'   $spearman    — Spearman correlation results (H_TRBV vs H_clonotype)
#'   $pearson     — Pearson correlation results
#'   $clonality_spearman — Spearman on clonality proxy vs true clonality
compare_diversity_metrics <- function(trbv_results,
                                      clonotype_results,
                                      group_vars = "sample_id") {

  joined <- inner_join(trbv_results, clonotype_results, by = group_vars)

  if (nrow(joined) < 3) {
    warning("Fewer than 3 matched groups; correlation not computed.")
    return(list(table = joined, spearman = NA, pearson = NA,
                clonality_spearman = NA))
  }

  # Shannon diversity correlation
  sp  <- cor.test(joined$H_TRBV, joined$H_clonotype,
                  method = "spearman", exact = FALSE)
  pe  <- cor.test(joined$H_TRBV, joined$H_clonotype,
                  method = "pearson")

  # Clonality correlation
  sp_clon <- cor.test(joined$TRBV_clonality_proxy, joined$true_clonality,
                      method = "spearman", exact = FALSE)

  list(
    table              = joined,
    spearman           = sp,
    pearson            = pe,
    clonality_spearman = sp_clon
  )
}


# =============================================================================
# SECTION 3: FAILURE MODE BENCHMARKING
# =============================================================================

#' Benchmark failure modes of the TRBV proxy
#'
#' For each group, reports:
#'  - How many clonotypes share each TRBV gene (within-TRBV clonotype count)
#'  - Cases where one expanded clone dominates within a TRBV family
#'  - Cases where TRBV diversity appears high but true clonotype diversity is low
#'  - Cases where TRBV diversity appears low but true clonotype diversity is high
#'
#' @param metadata    Data frame with columns: trbv_gene, clonotype_id, and
#'                    group_vars columns.
#' @param group_vars  Character vector of grouping columns.
#' @param min_cells   Minimum cells per group (default 30).
#' @return  List with:
#'   $within_trbv  — per-TRBV-gene clonotype statistics per group
#'   $failures     — joined diversity table flagged with failure mode labels
benchmark_failure_modes <- function(metadata,
                                    group_vars = "sample_id",
                                    min_cells  = 30) {

  req_cols <- c("trbv_gene", "clonotype_id")
  miss <- setdiff(req_cols, colnames(metadata))
  if (length(miss) > 0) stop("metadata missing columns: ", paste(miss, collapse = ", "))

  # 1. Within-TRBV clonotype diversity
  within_trbv <- metadata %>%
    filter(!is.na(trbv_gene), trbv_gene != "",
           !is.na(clonotype_id), clonotype_id != "") %>%
    group_by(across(all_of(c(group_vars, "trbv_gene")))) %>%
    summarise(
      n_cells_in_trbv     = n(),
      n_clonotypes_in_trbv = n_distinct(clonotype_id),
      # Shannon diversity of clonotypes within this TRBV gene
      H_within_trbv       = {
        cts <- table(clonotype_id)
        shannon_diversity(as.numeric(cts))
      },
      # Fraction of cells belonging to the dominant clone within this TRBV gene
      dominant_clone_fraction = {
        cts <- table(clonotype_id)
        max(cts) / sum(cts)
      },
      dominant_clone_id   = names(sort(table(clonotype_id), decreasing = TRUE))[1],
      .groups = "drop"
    )

  # 2. Flag samples where a single clone dominates a TRBV family (>50% of cells)
  expanded_within_trbv <- within_trbv %>%
    filter(dominant_clone_fraction > 0.5, n_cells_in_trbv >= 5) %>%
    arrange(desc(dominant_clone_fraction))

  # 3. Compute both diversity metrics and flag discordant cases
  trbv_div     <- compute_trbv_diversity(metadata, group_vars, min_cells)
  clono_div    <- compute_clonotype_diversity(metadata, group_vars, min_cells)
  joined       <- inner_join(trbv_div, clono_div, by = group_vars)

  if (nrow(joined) >= 2) {
    # Normalize H values to [0,1] for comparison
    joined <- joined %>%
      mutate(
        H_TRBV_z     = scale(H_TRBV)[,1],
        H_clono_z    = scale(H_clonotype)[,1],
        discordance  = H_TRBV_z - H_clono_z,
        failure_mode = case_when(
          discordance >  1 ~ "High TRBV diversity, Low true clonotype diversity",
          discordance < -1 ~ "Low TRBV diversity, High true clonotype diversity",
          TRUE             ~ "Concordant"
        )
      )
  }

  # 4. Compute how many clonotypes share each TRBV on average
  clono_per_trbv_summary <- within_trbv %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      mean_clonotypes_per_trbv = mean(n_clonotypes_in_trbv),
      max_clonotypes_per_trbv  = max(n_clonotypes_in_trbv),
      n_trbv_with_gt1_clono    = sum(n_clonotypes_in_trbv > 1),
      n_trbv_with_expanded     = sum(dominant_clone_fraction > 0.5 &
                                       n_cells_in_trbv >= 5),
      .groups = "drop"
    )

  list(
    within_trbv          = within_trbv,
    expanded_within_trbv = expanded_within_trbv,
    clono_per_trbv       = clono_per_trbv_summary,
    failures             = joined
  )
}


# =============================================================================
# SECTION 4: MAIN WRAPPER
# =============================================================================

#' Run complete TCR diversity analysis
#'
#' Computes TRBV usage diversity (NanoString-style proxy) and true clonotype
#' diversity, correlates them, and benchmarks failure modes.
#'
#' @param metadata     Data frame. Required columns: trbv_gene, clonotype_id, and
#'                     whatever columns are listed in group_vars.
#' @param group_vars   Character vector of grouping column names. Default "sample_id".
#' @param min_cells    Minimum cells per group (default 30).
#' @param downsample_n If > 0, repeat analysis on cells downsampled to this count.
#' @return  List with:
#'   $full          — results on all cells
#'   $downsampled   — results on downsampled cells (NULL if downsample_n = 0)
#'   $comparison    — output of compare_diversity_metrics()
#'   $benchmark     — output of benchmark_failure_modes()
run_tcr_diversity_analysis <- function(metadata,
                                       group_vars   = "sample_id",
                                       min_cells    = 30,
                                       downsample_n = 0) {

  message("=== TCR Diversity Analysis ===")
  message("Groups: ", paste(group_vars, collapse = " + "))
  message("Min cells per group: ", min_cells)

  # Full analysis
  message("\n[1/4] Computing TRBV usage diversity (NanoString-style proxy)...")
  trbv_full  <- compute_trbv_diversity(metadata, group_vars, min_cells)
  message("  -> ", nrow(trbv_full), " groups with TRBV data")

  message("[2/4] Computing true clonotype diversity...")
  clon_full  <- compute_clonotype_diversity(metadata, group_vars, min_cells)
  message("  -> ", nrow(clon_full), " groups with clonotype data")

  message("[3/4] Comparing TRBV proxy vs clonotype diversity...")
  comparison <- compare_diversity_metrics(trbv_full, clon_full, group_vars)
  n_matched  <- nrow(comparison$table)
  if (!is.null(comparison$spearman) && !all(is.na(comparison$spearman))) {
    rho <- round(comparison$spearman$estimate, 3)
    p   <- signif(comparison$spearman$p.value, 3)
    message("  -> Spearman rho (H_TRBV vs H_clonotype): ", rho,
            "  (p = ", p, ", n = ", n_matched, ")")
  }

  message("[4/4] Benchmarking failure modes...")
  benchmark  <- benchmark_failure_modes(metadata, group_vars, min_cells)

  # Optional downsampled analysis
  ds_results <- NULL
  if (downsample_n > 0) {
    message("\n[Sensitivity] Downsampling to ", downsample_n, " cells per group...")
    trbv_ds   <- compute_trbv_diversity(metadata, group_vars, min_cells, downsample_n)
    clon_ds   <- compute_clonotype_diversity(metadata, group_vars, min_cells, downsample_n)
    comp_ds   <- compare_diversity_metrics(trbv_ds, clon_ds, group_vars)
    ds_results <- list(trbv = trbv_ds, clonotype = clon_ds, comparison = comp_ds)
  }

  message("\nDone.")
  list(
    full       = list(trbv = trbv_full, clonotype = clon_full),
    comparison = comparison,
    benchmark  = benchmark,
    downsampled = ds_results
  )
}


# =============================================================================
# SECTION 5: DATA LOADING HELPERS
# =============================================================================

#' Load and combine VDJ data from multiple Cell Ranger per-sample outputs
#'
#' Reads filtered_contig_annotations.csv from each sample directory, keeps
#' only productive TRB contigs (one per cell — highest-UMI wins on ties),
#' and assigns gem-group-corrected barcodes to match cellranger aggr output.
#'
#' @param vdj_dir      Path to parent VDJ directory containing one sub-folder
#'                     per sample (e.g., "scVDJ/").
#' @param sample_names Character vector of sample names (sub-folder names),
#'                     in the same ORDER as they appear in the cellranger aggr
#'                     CSV (gem group 1 = sample_names[1], etc.).
#' @param file_pattern Glob pattern for the contig file within each sample dir.
#'                     Default = "*_filtered_contig_annotations.csv".
#' @return  Data frame with columns:
#'   barcode_aggr, sample_id, gem_group, chain, trbv_gene, clonotype_id,
#'   cdr3_aa, cdr3_nt, umis, is_cell, productive
load_vdj_data <- function(vdj_dir,
                           sample_names,
                           file_pattern = "*_filtered_contig_annotations.csv") {

  sample_files <- vapply(sample_names, function(s) {
    fls <- list.files(file.path(vdj_dir, s),
                      pattern = gsub("\\*", ".*", file_pattern),
                      full.names = TRUE)
    if (length(fls) == 0)
      stop("No contig file found for sample: ", s, " in ", file.path(vdj_dir, s))
    fls[1]
  }, character(1))

  purrr::imap_dfr(sample_files, function(fpath, sample_name) {
    gem_group <- which(sample_names == sample_name)

    df <- fread(fpath, data.table = FALSE, showProgress = FALSE) %>%
      # Standardise column names
      rename_with(~ tolower(gsub(" ", "_", .x)))

    # Keep only is_cell == TRUE rows (already filtered in filtered_ file, but safe)
    df <- df %>%
      filter(is_cell %in% c("true", TRUE, "True", 1))

    # Keep productive TRB chains only
    df_trb <- df %>%
      filter(chain == "TRB",
             productive %in% c("true", TRUE, "True", 1)) %>%
      # If a cell has multiple productive TRB contigs, keep highest-UMI one
      group_by(barcode) %>%
      slice_max(order_by = umis, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      transmute(
        # Replace the "-1" gem group suffix with the actual gem group index
        barcode_aggr = sub("-\\d+$", paste0("-", gem_group), barcode),
        sample_id    = sample_name,
        gem_group    = gem_group,
        chain        = chain,
        trbv_gene    = v_gene,
        trbd_gene    = if_else(d_gene == "" | is.na(d_gene), NA_character_, d_gene),
        trbj_gene    = j_gene,
        clonotype_id = paste0(sample_name, "_", raw_clonotype_id),
        cdr3_aa      = cdr3,
        cdr3_nt      = cdr3_nt,
        umis         = umis,
        raw_clonotype_id = raw_clonotype_id
      )

    # Also grab TRA contigs to build paired clonotype IDs (optional)
    df_tra <- df %>%
      filter(chain == "TRA",
             productive %in% c("true", TRUE, "True", 1)) %>%
      group_by(barcode) %>%
      slice_max(order_by = umis, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      transmute(
        barcode_aggr  = sub("-\\d+$", paste0("-", gem_group), barcode),
        trav_gene     = v_gene,
        cdr3_alpha_aa = cdr3
      )

    # Join TRA info where available
    df_trb <- left_join(df_trb, df_tra, by = "barcode_aggr")
    df_trb
  })
}

#' Read aggr barcodes and annotate with sample IDs
#'
#' @param aggr_dir     Path to cellranger aggr output directory containing
#'                     barcodes.tsv.
#' @param sample_names Character vector of sample names ordered by gem group.
#' @return  Data frame with columns: barcode_aggr, gem_group, sample_id.
load_aggr_barcodes <- function(aggr_dir, sample_names) {
  bc_file <- file.path(aggr_dir, "barcodes.tsv")
  if (!file.exists(bc_file)) stop("barcodes.tsv not found in: ", aggr_dir)

  barcodes <- fread(bc_file, header = FALSE, col.names = "barcode_aggr",
                    data.table = FALSE)
  barcodes <- barcodes %>%
    mutate(
      gem_group = as.integer(sub(".*-", "", barcode_aggr)),
      sample_id = sample_names[gem_group]
    )
  barcodes
}

#' Build a unified TCR metadata table from raw VDJ + aggr outputs
#'
#' Returns one row per cell barcode present in the aggr output that also has
#' a productive TRB contig. Cells without TCR data are dropped (they are not
#' T cells or lacked recovery).
#'
#' @param vdj_dir          Path to parent VDJ directory.
#' @param sample_names     Ordered character vector of sample names.
#' @param aggr_dir         Path to cellranger aggr output directory.
#' @param additional_meta  Optional data frame with extra per-cell metadata
#'                         (e.g., cell type annotations). Must contain
#'                         'barcode_aggr' as a join key.
#' @return  Metadata data frame ready for diversity analysis.
build_tcr_metadata <- function(vdj_dir,
                                sample_names,
                                aggr_dir          = NULL,
                                additional_meta   = NULL) {

  message("Loading VDJ data from ", length(sample_names), " samples...")
  vdj <- load_vdj_data(vdj_dir, sample_names)
  message("  -> ", nrow(vdj), " T cells with productive TRB contigs")

  if (!is.null(aggr_dir)) {
    message("Loading aggr barcodes from ", aggr_dir, "...")
    aggr_bc <- load_aggr_barcodes(aggr_dir, sample_names)
    message("  -> ", nrow(aggr_bc), " total cells in aggr")
    # Keep only cells that appear in the aggr output
    vdj <- inner_join(vdj, aggr_bc %>% select(barcode_aggr, gem_group, sample_id),
                      by = c("barcode_aggr", "gem_group", "sample_id"))
    message("  -> ", nrow(vdj), " T cells matched to aggr barcodes")
  }

  if (!is.null(additional_meta)) {
    if (!"barcode_aggr" %in% colnames(additional_meta)) {
      stop("additional_meta must contain column 'barcode_aggr'")
    }
    vdj <- left_join(vdj, additional_meta, by = "barcode_aggr")
  }

  vdj
}


# =============================================================================
# SECTION 6: PLOTTING FUNCTIONS
# =============================================================================

# Colour palette for TRBV genes (up to 40 distinct colours)
# Builds a qualitative palette by interpolating across multiple RColorBrewer sets.
.trbv_palette <- function(n) {
  if (n <= 8)  return(RColorBrewer::brewer.pal(max(n, 3), "Set2")[seq_len(n)])
  if (n <= 12) return(RColorBrewer::brewer.pal(12, "Set3")[seq_len(n)])
  # For larger n, interpolate across Paired + Set1 + Dark2
  base_cols <- c(
    RColorBrewer::brewer.pal(12, "Paired"),
    RColorBrewer::brewer.pal(8,  "Dark2"),
    RColorBrewer::brewer.pal(9,  "Set1")
  )
  grDevices::colorRampPalette(base_cols)(n)
}

#' Stacked barplot of TRBV gene usage per group
#'
#' @param metadata   Data frame with columns: trbv_gene, and group_var column.
#' @param group_var  Single column name to use on the x-axis (e.g., "sample_id").
#' @param top_n_genes  Show only top N TRBV genes; remainder grouped as "Other".
#' @return  ggplot2 object.
plot_trbv_stacked_bar <- function(metadata,
                                   group_var   = "sample_id",
                                   top_n_genes = 20) {

  # Identify top TRBV genes by total usage across all groups
  top_genes <- metadata %>%
    filter(!is.na(trbv_gene), trbv_gene != "") %>%
    count(trbv_gene, sort = TRUE) %>%
    slice_head(n = top_n_genes) %>%
    pull(trbv_gene)

  plot_df <- metadata %>%
    filter(!is.na(trbv_gene), trbv_gene != "") %>%
    mutate(trbv_plot = if_else(trbv_gene %in% top_genes, trbv_gene, "Other")) %>%
    count(.data[[group_var]], trbv_plot) %>%
    group_by(.data[[group_var]]) %>%
    mutate(freq = n / sum(n)) %>%
    ungroup() %>%
    mutate(trbv_plot = factor(trbv_plot,
                              levels = c(sort(top_genes), "Other")))

  n_colors <- length(unique(plot_df$trbv_plot))
  colors   <- .trbv_palette(n_colors - 1)  # -1 for "Other"
  colors   <- c(colors, "grey70")
  names(colors) <- levels(plot_df$trbv_plot)

  ggplot(plot_df, aes(x = .data[[group_var]], y = freq, fill = trbv_plot)) +
    geom_bar(stat = "identity", width = 0.8, color = "white", linewidth = 0.2) +
    scale_fill_manual(values = colors, name = "TRBV gene") +
    scale_y_continuous(labels = percent_format()) +
    labs(
      title    = "TRBV Gene Usage per Sample",
      subtitle = paste0("NanoString-style TRBV usage diversity proxy | Top ",
                        top_n_genes, " genes shown"),
      x        = group_var,
      y        = "Proportion of T cells"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      legend.key.size = unit(0.4, "cm"),
      legend.text  = element_text(size = 8)
    )
}

#' Clonotype rank-abundance curves (Zipf plot) per group
#'
#' @param metadata  Data frame with columns: clonotype_id, and group_var column.
#' @param group_var Single column name (e.g., "sample_id").
#' @param top_n     Number of top clonotypes to label.
#' @return  ggplot2 object.
plot_clonotype_rank_abundance <- function(metadata,
                                          group_var = "sample_id",
                                          top_n     = 5) {

  rank_df <- metadata %>%
    filter(!is.na(clonotype_id), clonotype_id != "") %>%
    group_by(.data[[group_var]], clonotype_id) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(.data[[group_var]]) %>%
    arrange(desc(n)) %>%
    mutate(
      rank = row_number(),
      freq = n / sum(n)
    ) %>%
    ungroup()

  ggplot(rank_df, aes(x = rank, y = freq,
                       color = .data[[group_var]],
                       group = .data[[group_var]])) +
    geom_line(alpha = 0.8, linewidth = 0.8) +
    scale_x_log10(labels = label_number()) +
    scale_y_log10(labels = percent_format()) +
    labs(
      title    = "Clonotype Rank-Abundance Curves",
      subtitle = "Rank 1 = most abundant clonotype; steeper curve = more clonal",
      x        = "Clonotype rank (log scale)",
      y        = "Proportion of T cells (log scale)",
      color    = group_var
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")
}

#' Scatterplot of TRBV usage diversity vs true clonotype diversity
#'
#' @param diversity_table  Output of compare_diversity_metrics()$table.
#' @param group_var        Column used to colour/label points.
#' @param label_col        Column to use for point labels (default = group_var).
#' @return  ggplot2 object.
plot_h_trbv_vs_h_clonotype <- function(diversity_table,
                                        group_var  = "sample_id",
                                        label_col  = NULL) {

  if (is.null(label_col)) label_col <- group_var

  # Correlation annotation
  sp <- tryCatch(
    cor.test(diversity_table$H_TRBV, diversity_table$H_clonotype,
             method = "spearman", exact = FALSE),
    error = function(e) NULL
  )
  subtitle <- if (!is.null(sp)) {
    paste0("Spearman rho = ", round(sp$estimate, 3),
           "  (p = ", signif(sp$p.value, 3), ")")
  } else "Insufficient data for correlation"

  ggplot(diversity_table,
         aes(x = H_TRBV, y = H_clonotype, color = .data[[group_var]])) +
    geom_point(size = 3, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, color = "grey40",
                linetype = "dashed", linewidth = 0.8) +
    geom_text_repel(aes(label = .data[[label_col]]),
                    size = 3, max.overlaps = 15) +
    geom_abline(slope = 1, intercept = 0,
                color = "red", linetype = "dotted", linewidth = 0.6) +
    labs(
      title    = "TRBV Usage Diversity vs True Clonotype Diversity",
      subtitle = subtitle,
      x        = expression(H[TRBV]~"(TRBV usage diversity proxy)"),
      y        = expression(H[clonotype]~"(true CDR3-level diversity)"),
      caption  = paste("Red dotted line: 1:1 reference.",
                       "TRBV proxy cannot resolve CDR3 clonotypes.")
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")
}

#' Scatterplot of TRBV clonality proxy vs true clonality
#'
#' @param diversity_table  Output of compare_diversity_metrics()$table.
#' @param group_var        Column used to colour/label points.
#' @return  ggplot2 object.
plot_clonality_comparison <- function(diversity_table,
                                       group_var = "sample_id") {

  sp <- tryCatch(
    cor.test(diversity_table$TRBV_clonality_proxy, diversity_table$true_clonality,
             method = "spearman", exact = FALSE),
    error = function(e) NULL
  )
  subtitle <- if (!is.null(sp)) {
    paste0("Spearman rho = ", round(sp$estimate, 3),
           "  (p = ", signif(sp$p.value, 3), ")")
  } else "Insufficient data for correlation"

  ggplot(diversity_table,
         aes(x = TRBV_clonality_proxy, y = true_clonality,
             color = .data[[group_var]])) +
    geom_point(size = 3, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, color = "grey40",
                linetype = "dashed", linewidth = 0.8) +
    geom_text_repel(aes(label = .data[[group_var]]),
                    size = 3, max.overlaps = 15) +
    geom_abline(slope = 1, intercept = 0,
                color = "red", linetype = "dotted", linewidth = 0.6) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title    = "TRBV Clonality Proxy vs True Clonality",
      subtitle = subtitle,
      x        = "TRBV clonality proxy (1 - H_norm_TRBV)",
      y        = "True clonality (1 - H_norm_clonotype)"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")
}

#' Heatmap of dominant TRBV genes by sample (proportion matrix)
#'
#' @param metadata    Data frame with trbv_gene and group_var columns.
#' @param group_var   Column for columns of the heatmap (e.g., "sample_id").
#' @param top_n_genes Number of top TRBV genes to display (default 25).
#' @return  Invisible NULL (pheatmap renders directly).
plot_trbv_heatmap <- function(metadata,
                               group_var   = "sample_id",
                               top_n_genes = 25) {

  top_genes <- metadata %>%
    filter(!is.na(trbv_gene), trbv_gene != "") %>%
    count(trbv_gene, sort = TRUE) %>%
    slice_head(n = top_n_genes) %>%
    pull(trbv_gene)

  heat_mat <- metadata %>%
    filter(trbv_gene %in% top_genes) %>%
    count(.data[[group_var]], trbv_gene) %>%
    group_by(.data[[group_var]]) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    select(-n) %>%
    pivot_wider(names_from = all_of(group_var), values_from = prop,
                values_fill = 0) %>%
    tibble::column_to_rownames("trbv_gene") %>%
    as.matrix()

  pheatmap::pheatmap(
    heat_mat,
    scale        = "none",
    color        = colorRampPalette(c("white", "#2166AC", "#053061"))(100),
    main         = paste("TRBV Gene Usage Heatmap\n(proportion of cells per sample)"),
    fontsize_row = 8,
    fontsize_col = 10,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    border_color = NA
  )
  invisible(NULL)
}

#' Plot within-TRBV clonotype diversity (failure mode benchmark)
#'
#' Shows how many clonotypes reside within each TRBV gene family, highlighting
#' genes where a single clone dominates.
#'
#' @param benchmark_results  Output of benchmark_failure_modes().
#' @param group_var          Column used to facet (e.g., "sample_id").
#' @param top_n_trbv         Number of top TRBV genes to show.
#' @return  ggplot2 object.
plot_within_trbv_diversity <- function(benchmark_results,
                                        group_var  = "sample_id",
                                        top_n_trbv = 15) {

  df <- benchmark_results$within_trbv

  top_trbv <- df %>%
    group_by(trbv_gene) %>%
    summarise(total = sum(n_cells_in_trbv), .groups = "drop") %>%
    slice_max(total, n = top_n_trbv) %>%
    pull(trbv_gene)

  df_plot <- df %>%
    filter(trbv_gene %in% top_trbv) %>%
    mutate(
      trbv_gene = factor(trbv_gene, levels = rev(top_trbv)),
      expansion = dominant_clone_fraction > 0.5
    )

  ggplot(df_plot,
         aes(x = n_clonotypes_in_trbv, y = trbv_gene,
             fill = dominant_clone_fraction)) +
    geom_bar(stat = "identity") +
    facet_wrap(vars(.data[[group_var]]), scales = "free_x") +
    scale_fill_gradient2(
      low      = "#4DAF4A",
      mid      = "#FF7F00",
      high     = "#E41A1C",
      midpoint = 0.5,
      limits   = c(0, 1),
      name     = "Dominant\nclone\nfraction"
    ) +
    labs(
      title    = "Within-TRBV Clonotype Diversity (Failure Mode Analysis)",
      subtitle = "Red fill: single clone >50% of cells within that TRBV gene family",
      x        = "Number of unique clonotypes within TRBV gene",
      y        = "TRBV gene"
    ) +
    theme_bw(base_size = 10) +
    theme(strip.text = element_text(size = 8))
}

#' Bland-Altman style plot: difference vs mean of TRBV and clonotype diversity
#'
#' Visualises systematic bias and limits of agreement between the proxy and
#' the true diversity metric.
#'
#' @param diversity_table  Output of compare_diversity_metrics()$table.
#' @param group_var        Column for point labels.
#' @return  ggplot2 object.
plot_bland_altman <- function(diversity_table, group_var = "sample_id") {

  ba_df <- diversity_table %>%
    mutate(
      mean_H = (H_TRBV + H_clonotype) / 2,
      diff_H = H_TRBV - H_clonotype
    )

  mean_diff <- mean(ba_df$diff_H, na.rm = TRUE)
  sd_diff   <- sd(ba_df$diff_H, na.rm = TRUE)
  loa_upper <- mean_diff + 1.96 * sd_diff
  loa_lower <- mean_diff - 1.96 * sd_diff

  ggplot(ba_df, aes(x = mean_H, y = diff_H, color = .data[[group_var]])) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = .data[[group_var]]), size = 3) +
    geom_hline(yintercept = mean_diff, color = "blue",  linetype = "solid") +
    geom_hline(yintercept = loa_upper, color = "red",   linetype = "dashed") +
    geom_hline(yintercept = loa_lower, color = "red",   linetype = "dashed") +
    geom_hline(yintercept = 0,         color = "grey50", linetype = "dotted") +
    annotate("text", x = max(ba_df$mean_H, na.rm = TRUE),
             y = mean_diff + 0.05, label = paste("Mean diff:", round(mean_diff, 3)),
             hjust = 1, size = 3, color = "blue") +
    annotate("text", x = max(ba_df$mean_H, na.rm = TRUE),
             y = loa_upper + 0.05, label = paste("+1.96 SD:", round(loa_upper, 3)),
             hjust = 1, size = 3, color = "red") +
    annotate("text", x = max(ba_df$mean_H, na.rm = TRUE),
             y = loa_lower - 0.05, label = paste("-1.96 SD:", round(loa_lower, 3)),
             hjust = 1, size = 3, color = "red") +
    labs(
      title    = "Bland-Altman Plot: H_TRBV vs H_clonotype",
      subtitle = "Bias and limits of agreement between TRBV proxy and true diversity",
      x        = "Mean of H_TRBV and H_clonotype",
      y        = "H_TRBV - H_clonotype (difference)"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")
}


# =============================================================================
# SECTION 7: SEURAT INTEGRATION HELPERS
# =============================================================================

#' Extract metadata from a Seurat object and attach TCR data
#'
#' Merges Seurat cell-level metadata with VDJ data by cell barcode. Optionally
#' filters to T cells using canonical marker gene expression thresholds.
#'
#' @param seurat_obj       A Seurat object (Seurat v4 or v5).
#' @param vdj_data         Data frame from load_vdj_data() with barcode_aggr column.
#' @param tcell_annotation Column name in Seurat metadata containing cell type
#'                         annotations. Cells whose value contains "T" (case-
#'                         insensitive) are kept. Set to NULL to skip this filter.
#' @param cd3_threshold    If > 0, require CD3D + CD3E + CD3G normalised
#'                         expression sum >= this threshold (in addition to or
#'                         instead of annotation filter).
#' @return  Data frame with merged Seurat metadata + TCR columns.
extract_seurat_metadata <- function(seurat_obj,
                                     vdj_data,
                                     tcell_annotation = "cell_type",
                                     cd3_threshold    = 0) {

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required. Install with: install.packages('Seurat')")
  }

  # Pull cell metadata
  meta <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("barcode_aggr")

  # Optional T cell annotation filter
  if (!is.null(tcell_annotation) && tcell_annotation %in% colnames(meta)) {
    meta <- meta %>%
      filter(grepl("T cell|T-cell|^T$|CD4|CD8|Treg|NK-T|NKT",
                   .data[[tcell_annotation]], ignore.case = TRUE))
    message("  -> ", nrow(meta), " annotated T cells retained")
  }

  # Optional CD3 expression filter
  if (cd3_threshold > 0) {
    cd3_genes <- intersect(c("CD3D", "CD3E", "CD3G"),
                           rownames(seurat_obj))
    if (length(cd3_genes) == 0) {
      warning("No CD3D/CD3E/CD3G genes found in Seurat object; skipping CD3 filter.")
    } else {
      cd3_exp <- Matrix::colSums(
        Seurat::GetAssayData(seurat_obj, layer = "data")[cd3_genes, , drop = FALSE]
      )
      cd3_cells <- names(cd3_exp[cd3_exp >= cd3_threshold])
      meta <- meta %>% filter(barcode_aggr %in% cd3_cells)
      message("  -> ", nrow(meta), " CD3+ cells retained")
    }
  }

  # Merge with VDJ
  merged <- left_join(meta, vdj_data, by = "barcode_aggr")
  message("  -> ", sum(!is.na(merged$trbv_gene)), " cells with TRBV assignment")
  merged
}


# =============================================================================
# SECTION 8: OUTPUT FORMATTING
# =============================================================================

#' Build the final summary table combining TRBV and clonotype diversity results
#'
#' @param comparison  Output of compare_diversity_metrics().
#' @return  Formatted tibble with all output columns.
build_summary_table <- function(comparison) {
  comparison$table %>%
    select(
      any_of(c("sample_id", "subset", "condition", "patient_id")),
      n_cells_trbv,
      n_trbv_detected,
      H_TRBV,
      H_norm_TRBV,
      TRBV_clonality_proxy,
      n_clonotypes,
      H_clonotype,
      H_norm_clonotype,
      true_clonality,
      top_clone_fraction
    ) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4)))
}

#' Print a formatted correlation summary
#'
#' @param comparison  Output of compare_diversity_metrics().
print_correlation_summary <- function(comparison) {
  cat("\n=== Diversity Metric Correlations ===\n")
  cat("(TRBV usage diversity proxy vs true CDR3 clonotype diversity)\n\n")

  if (!is.null(comparison$spearman) && !all(is.na(comparison$spearman))) {
    sp  <- comparison$spearman
    pe  <- comparison$pearson
    sp2 <- comparison$clonality_spearman

    cat(sprintf("Shannon diversity (H_TRBV vs H_clonotype):\n"))
    cat(sprintf("  Spearman rho = %.3f  (p = %.4f)\n",
                sp$estimate, sp$p.value))
    cat(sprintf("  Pearson  r   = %.3f  (p = %.4f)\n\n",
                pe$estimate, pe$p.value))
    cat(sprintf("Clonality (TRBV proxy vs true clonality):\n"))
    cat(sprintf("  Spearman rho = %.3f  (p = %.4f)\n\n",
                sp2$estimate, sp2$p.value))
  } else {
    cat("Insufficient data for correlation analysis.\n")
  }
  invisible(NULL)
}


# =============================================================================
# SECTION 9: T CELL SUBSET ANNOTATION FROM scRNA-seq EXPRESSION
# =============================================================================
#
# PURPOSE:
#   Classify cells in a cellranger aggr output into CD4, CD8, Treg, and
#   double-negative (DN) T cell compartments using marker gene expression.
#   This enables compartment-stratified diversity analysis — the key use case
#   for NanoString-style TRBV diversity: estimating clonal structure within
#   a defined T cell subset (e.g. sorted CD8 T cells from an nCounter run).
#
# MARKER LOGIC (applied in hierarchical order):
#   CD8   :  CD8A > 0  OR  CD8B > 0
#   Treg  :  FOXP3 > 0  AND CD4 > 0  AND NOT CD8
#   CD4   :  CD4 > 0  AND NOT CD8  (includes Treg unless split = TRUE)
#   DN    :  CD3D > 0  OR CD3E > 0  (no CD4 or CD8 signal)
#   Non-T :  no CD3 signal
#
# NOTE: This is a conservative marker-gene-based classification, not a full
# clustering pipeline. It works well for clean CD4/CD8 separation but will
# miss cells with low expression or ambient RNA contamination. Use a proper
# Seurat clustering + marker annotation workflow for publication-quality
# subset assignments.
# =============================================================================

#' Annotate T cell subsets from a cellranger aggr sparse expression matrix
#'
#' Reads the aggr matrix.mtx file and classifies each barcode as CD8, CD4,
#' Treg, DN (double-negative), or Non-T using marker gene expression.
#' Uses Matrix::readMM for correct, efficient sparse matrix loading.
#'
#' @param aggr_dir    Path to cellranger aggr output directory (contains
#'                    matrix.mtx, features.tsv, barcodes.tsv).
#' @param count_threshold  Minimum raw UMI count to call a marker positive.
#'                         Default = 1 (any detected count).
#' @param split_treg  If TRUE, Tregs are labeled "Treg" separately from "CD4".
#'                    If FALSE (default), Tregs are grouped with "CD4".
#' @return  data.table with columns: barcode_aggr, subset, cd4_count, cd8a_count,
#'          cd8b_count, cd3d_count, cd3e_count, foxp3_count.
annotate_tcell_subsets_from_aggr <- function(aggr_dir,
                                              count_threshold = 1,
                                              split_treg      = FALSE) {

  mtx_file  <- file.path(aggr_dir, "matrix.mtx")
  feat_file <- file.path(aggr_dir, "features.tsv")
  bc_file   <- file.path(aggr_dir, "barcodes.tsv")

  for (f in c(mtx_file, feat_file, bc_file)) {
    if (!file.exists(f)) stop("File not found: ", f)
  }

  message("  Reading expression matrix (Matrix::readMM)...")
  mat <- Matrix::readMM(mtx_file)   # features x cells sparse matrix

  # Find row indices of marker genes in the features file
  feat <- fread(feat_file, header = FALSE,
                col.names = c("ensembl", "symbol", "type"), showProgress = FALSE)
  feat[, row_idx := .I]

  markers <- c("CD4", "CD8A", "CD8B", "CD3D", "CD3E", "FOXP3")
  m_tbl   <- feat[symbol %in% markers, .(symbol, row_idx)]
  missing <- setdiff(markers, m_tbl$symbol)
  if (length(missing) > 0) {
    warning("Marker genes not found in features.tsv: ", paste(missing, collapse = ", "))
  }

  # Subset matrix to marker gene rows only
  sub <- mat[m_tbl$row_idx, , drop = FALSE]
  rownames(sub) <- m_tbl$symbol

  # Load barcodes
  bc <- fread(bc_file, header = FALSE, col.names = "barcode_aggr",
              showProgress = FALSE)

  # Build per-cell count data.table
  get_row <- function(gene) {
    if (gene %in% rownames(sub)) as.integer(sub[gene, ]) else integer(ncol(sub))
  }
  cell_dt <- data.table(
    barcode_aggr = bc$barcode_aggr,
    cd4_count    = get_row("CD4"),
    cd8a_count   = get_row("CD8A"),
    cd8b_count   = get_row("CD8B"),
    cd3d_count   = get_row("CD3D"),
    cd3e_count   = get_row("CD3E"),
    foxp3_count  = get_row("FOXP3")
  )

  thr <- count_threshold

  # Hierarchical classification
  cell_dt[, subset := fcase(
    (cd8a_count >= thr | cd8b_count >= thr),                         "CD8",
    (cd4_count  >= thr & (cd8a_count < thr & cd8b_count < thr) &
       foxp3_count >= thr),                                           "Treg",
    (cd4_count  >= thr & (cd8a_count < thr & cd8b_count < thr)),     "CD4",
    (cd3d_count >= thr | cd3e_count >= thr),                         "DN_T",
    default = "Non_T"
  )]

  # Optionally collapse Treg into CD4
  if (!split_treg) {
    cell_dt[subset == "Treg", subset := "CD4"]
  }

  # Summary
  summary_tbl <- cell_dt[, .N, by = subset][order(-N)]
  message("  T cell subset annotation complete:")
  for (i in seq_len(nrow(summary_tbl))) {
    message(sprintf("    %-8s: %d cells (%.1f%%)",
                    summary_tbl$subset[i], summary_tbl$N[i],
                    100 * summary_tbl$N[i] / nrow(cell_dt)))
  }

  cell_dt[]
}

#' Join T cell subset annotations onto a VDJ metadata table
#'
#' @param vdj_metadata   Data frame from build_tcr_metadata(), with barcode_aggr.
#' @param subset_annot   Data table from annotate_tcell_subsets_from_aggr().
#' @param keep_non_t     If FALSE (default), drop cells classified as Non_T.
#' @return  Merged metadata with a 'subset' column added.
add_subset_annotation <- function(vdj_metadata,
                                   subset_annot,
                                   keep_non_t = FALSE) {

  merged <- left_join(
    vdj_metadata,
    as.data.frame(subset_annot)[, c("barcode_aggr", "subset",
                                     "cd4_count", "cd8a_count",
                                     "cd8b_count", "foxp3_count")],
    by = "barcode_aggr"
  )

  # Cells in VDJ data but not in subset_annot get NA subset — label as "Unknown"
  merged$subset[is.na(merged$subset)] <- "Unknown"

  if (!keep_non_t) {
    n_before <- nrow(merged)
    merged   <- merged[merged$subset != "Non_T", ]
    n_drop   <- n_before - nrow(merged)
    if (n_drop > 0) message("  Dropped ", n_drop, " Non_T cells from VDJ metadata")
  }

  tab <- table(merged$subset)
  message("  VDJ-matched cells per subset:")
  for (nm in names(tab)) message(sprintf("    %-8s: %d", nm, tab[nm]))

  merged
}


# =============================================================================
# SECTION 10: MULTI-SEGMENT DIVERSITY (TRBV + TRBD + TRBJ + COMBINATIONS)
# =============================================================================
#
# RATIONALE:
#   The NanoString TRBV proxy only uses V-gene usage. Here we extend the analysis
#   to TRBD (D gene), TRBJ (J gene), and joint VDJ/VJ gene-combination usage as
#   richer proxies. Each combination narrows the resolution gap toward true CDR3
#   clonotype diversity, though none can fully resolve it.
#
#   NOTE ON TRBD:
#   In the J1568 dataset, ~83% of TRB contigs have no assigned D gene.
#   D-segment assignment is inherently unreliable due to the short segment length
#   and extensive junctional trimming. TRBD-based diversity should be interpreted
#   with extreme caution; VJ and VDJ combo metrics fall back to VJ when D is absent.
#
# KEY METRICS (in order of increasing resolution):
#   TRBV only         (V alone,  44 genes)
#   TRBJ only         (J alone,  13 genes)
#   TRBD only         (D alone,   2 genes — low information, high missingness)
#   VJ  combination   (V×J,      440 unique in J1568)
#   VDJ combination   (V×D×J,    687 unique, D-assigned cells only)
#   True clonotype    (CDR3,    2866 unique — ground truth)
# =============================================================================

#' Compute Shannon diversity for an arbitrary gene-usage column
#'
#' Generalised version of compute_trbv_diversity() that works on any named
#' gene column (trbv_gene, trbj_gene, trbd_gene, vj_combo, vdj_combo, ...).
#'
#' @param metadata   Data frame with the gene column and group_vars columns.
#' @param gene_col   Name of the gene/combo column to compute diversity on.
#' @param group_vars Character vector of grouping column names.
#' @param min_cells  Minimum cells per group (applied AFTER NA removal).
#' @param downsample_n  If > 0, downsample before computing.
#' @param metric_label  Short label used in the output column names
#'                      (default = gene_col).
#' @return  Data frame with one row per group and diversity columns prefixed
#'          by metric_label (e.g., H_trbj_gene, H_norm_trbj_gene, ...).
compute_gene_segment_diversity <- function(metadata,
                                            gene_col,
                                            group_vars    = "sample_id",
                                            min_cells     = 30,
                                            downsample_n  = 0,
                                            metric_label  = gene_col) {

  if (!gene_col %in% colnames(metadata)) {
    stop("Column '", gene_col, "' not found in metadata.")
  }
  if (downsample_n > 0) {
    metadata <- downsample_cells(metadata, group_vars, downsample_n)
  }

  result <- metadata %>%
    filter(!is.na(.data[[gene_col]]), .data[[gene_col]] != "",
           .data[[gene_col]] != "None") %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      n_cells          = n(),
      n_detected       = n_distinct(.data[[gene_col]]),
      H                = { cts <- table(.data[[gene_col]]); shannon_diversity(as.numeric(cts)) },
      H_norm           = { cts <- table(.data[[gene_col]]); normalized_shannon(as.numeric(cts)) },
      clonality_proxy  = { cts <- table(.data[[gene_col]]); clonality_from_shannon(as.numeric(cts)) },
      simpson          = { cts <- table(.data[[gene_col]]); simpson_diversity(as.numeric(cts)) },
      top_gene_fraction = { cts <- table(.data[[gene_col]]); max(cts) / sum(cts) },
      .groups = "drop"
    ) %>%
    filter(n_cells >= min_cells)

  # Rename generic columns to metric-specific names
  result <- result %>%
    rename(
      !!paste0("n_cells_",     metric_label) := n_cells,
      !!paste0("n_detected_",  metric_label) := n_detected,
      !!paste0("H_",           metric_label) := H,
      !!paste0("H_norm_",      metric_label) := H_norm,
      !!paste0("clonality_",   metric_label) := clonality_proxy,
      !!paste0("simpson_",     metric_label) := simpson,
      !!paste0("top_frac_",    metric_label) := top_gene_fraction
    )
  result
}

#' Build VJ and VDJ gene combination columns in metadata
#'
#' Adds two new columns:
#'   vj_combo  — paste(trbv_gene, trbj_gene, sep = ":")
#'   vdj_combo — paste(trbv_gene, trbd_gene, trbj_gene, sep = ":"), NA if D missing
#'
#' @param metadata  Data frame with trbv_gene, trbd_gene, trbj_gene columns.
#' @return  metadata with vj_combo and vdj_combo columns added.
add_vdj_combos <- function(metadata) {
  if (!all(c("trbv_gene", "trbj_gene") %in% colnames(metadata))) {
    stop("metadata must contain trbv_gene and trbj_gene (run load_vdj_data first).")
  }
  metadata <- metadata %>%
    mutate(
      # VJ combo: always defined when both V and J are assigned
      vj_combo = if_else(
        !is.na(trbv_gene) & !is.na(trbj_gene) & trbv_gene != "" & trbj_gene != "",
        paste(trbv_gene, trbj_gene, sep = ":"),
        NA_character_
      ),
      # VDJ combo: only defined when D is also assigned (only ~17% of cells in J1568)
      vdj_combo = if_else(
        !is.na(trbv_gene) & !is.na(trbd_gene) & !is.na(trbj_gene) &
          trbv_gene != "" & trbd_gene != "" & trbj_gene != "",
        paste(trbv_gene, trbd_gene, trbj_gene, sep = ":"),
        NA_character_
      )
    )
  metadata
}

#' Run all gene-segment diversity metrics and compare each to true clonotype diversity
#'
#' Computes Shannon diversity for TRBV, TRBJ, TRBD, VJ combos, and VDJ combos,
#' then correlates each with true clonotype diversity.
#'
#' @param metadata    Data frame from build_tcr_metadata() or add_vdj_combos().
#'                    Must contain: trbv_gene, trbj_gene, trbd_gene, clonotype_id,
#'                    vj_combo, vdj_combo (run add_vdj_combos() first).
#' @param group_vars  Character vector of grouping columns.
#' @param min_cells   Minimum cells per group (default 30).
#' @param downsample_n  Optional downsampling target.
#' @return  List with:
#'   $metrics      — wide table, one row per group, all diversity columns
#'   $correlations — data frame of Spearman rho, Pearson r, p-value per metric
#'   $clonotype    — clonotype diversity table (for reference)
run_multisegment_analysis <- function(metadata,
                                       group_vars   = "sample_id",
                                       min_cells    = 30,
                                       downsample_n = 0) {

  # Add VJ/VDJ combo columns if not already present
  if (!"vj_combo" %in% colnames(metadata)) {
    metadata <- add_vdj_combos(metadata)
  }

  message("=== Multi-Segment Diversity Analysis ===")
  n_d_assigned <- sum(!is.na(metadata$trbd_gene))
  n_total      <- nrow(metadata)
  message(sprintf("  TRBD assigned in %.0f%% of cells (%d / %d)",
                  100 * n_d_assigned / n_total, n_d_assigned, n_total))

  # Define the metrics to run (col_name, label displayed in output)
  segments <- list(
    list(col = "trbv_gene",  label = "TRBV"),
    list(col = "trbj_gene",  label = "TRBJ"),
    list(col = "trbd_gene",  label = "TRBD"),
    list(col = "vj_combo",   label = "VJ"),
    list(col = "vdj_combo",  label = "VDJ")
  )

  # Compute clonotype diversity (the reference)
  clono <- compute_clonotype_diversity(metadata, group_vars, min_cells, downsample_n)

  # Compute each segment metric and join sequentially
  seg_tables <- lapply(segments, function(seg) {
    # Skip if column missing from data
    if (!seg$col %in% colnames(metadata)) {
      message("  Skipping ", seg$label, " (column not found)")
      return(NULL)
    }
    message("  Computing ", seg$label, " diversity...")
    compute_gene_segment_diversity(
      metadata, seg$col, group_vars, min_cells, downsample_n,
      metric_label = seg$label
    )
  })
  # Drop NULLs
  seg_tables <- Filter(Negate(is.null), seg_tables)

  # Join all segment tables together by group_vars
  metrics_wide <- Reduce(function(a, b) left_join(a, b, by = group_vars),
                          seg_tables)
  # Join clonotype diversity
  metrics_wide <- left_join(metrics_wide, clono, by = group_vars)

  # Compute correlations of each metric's H with H_clonotype
  h_cols <- paste0("H_", vapply(segments, `[[`, character(1), "label"))
  h_cols <- intersect(h_cols, colnames(metrics_wide))  # keep only present ones

  correlations <- purrr::map_dfr(h_cols, function(hcol) {
    label <- sub("^H_", "", hcol)
    x <- metrics_wide[[hcol]]
    y <- metrics_wide$H_clonotype

    # Require at least 3 complete pairs
    complete <- !is.na(x) & !is.na(y)
    if (sum(complete) < 3) {
      return(data.frame(
        metric = label, n_groups = sum(complete),
        spearman_rho = NA_real_, spearman_p = NA_real_,
        pearson_r = NA_real_,    pearson_p  = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    sp <- cor.test(x[complete], y[complete], method = "spearman", exact = FALSE)
    pe <- cor.test(x[complete], y[complete], method = "pearson")
    data.frame(
      metric       = label,
      n_groups     = sum(complete),
      spearman_rho = round(sp$estimate, 4),
      spearman_p   = round(sp$p.value,  4),
      pearson_r    = round(pe$estimate, 4),
      pearson_p    = round(pe$p.value,  4),
      stringsAsFactors = FALSE
    )
  })

  # Also compute clonality correlations
  clon_cols <- paste0("clonality_", vapply(segments, `[[`, character(1), "label"))
  clon_cols <- intersect(clon_cols, colnames(metrics_wide))

  clon_correlations <- purrr::map_dfr(clon_cols, function(ccol) {
    label <- sub("^clonality_", "", ccol)
    x <- metrics_wide[[ccol]]
    y <- metrics_wide$true_clonality
    complete <- !is.na(x) & !is.na(y)
    if (sum(complete) < 3) return(NULL)
    sp <- cor.test(x[complete], y[complete], method = "spearman", exact = FALSE)
    data.frame(
      metric = label, n_groups = sum(complete),
      spearman_rho_clonality = round(sp$estimate, 4),
      spearman_p_clonality   = round(sp$p.value,  4),
      stringsAsFactors = FALSE
    )
  })

  correlations <- left_join(correlations, clon_correlations, by = c("metric", "n_groups"))
  correlations <- correlations %>% arrange(desc(abs(spearman_rho)))

  message("\nCorrelation with true clonotype diversity (Spearman rho, H metric):")
  print(correlations %>% select(metric, n_groups, spearman_rho, spearman_p))

  list(
    metrics      = metrics_wide,
    correlations = correlations,
    clonotype    = clono
  )
}

#' Faceted scatterplot: each gene-segment diversity metric vs H_clonotype
#'
#' @param ms_results  Output of run_multisegment_analysis().
#' @param group_var   Column used for point labels.
#' @return  ggplot2 object.
plot_metric_correlation_panel <- function(ms_results, group_var = "sample_id") {

  wide <- ms_results$metrics

  # Collect H_ columns for each segment metric
  h_cols <- grep("^H_(TRBV|TRBJ|TRBD|VJ|VDJ)$", colnames(wide), value = TRUE)

  # Build long-format data frame
  plot_df <- purrr::map_dfr(h_cols, function(hcol) {
    label <- sub("^H_", "", hcol)
    wide %>%
      select(all_of(c(group_vars_from_wide(wide, group_var), hcol, "H_clonotype"))) %>%
      rename(H_metric = all_of(hcol)) %>%
      mutate(metric = label)
  })

  # Correlation labels per facet
  rho_labels <- ms_results$correlations %>%
    transmute(
      metric,
      label = paste0("rho = ", sprintf("%.2f", spearman_rho),
                     "\np = ",  sprintf("%.3f", spearman_p))
    )

  # Axis range for 1:1 reference line
  all_h <- c(plot_df$H_metric, plot_df$H_clonotype)
  h_range <- range(all_h, na.rm = TRUE)

  ggplot(plot_df, aes(x = H_metric, y = H_clonotype)) +
    geom_point(aes(color = .data[[group_var]]), size = 2.5, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, color = "grey40",
                linetype = "dashed", linewidth = 0.7) +
    geom_text_repel(aes(label = .data[[group_var]]),
                    size = 2.5, max.overlaps = 10) +
    geom_abline(slope = 1, intercept = 0,
                color = "red", linetype = "dotted", linewidth = 0.5) +
    geom_text(data = rho_labels,
              aes(label = label),
              x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2,
              size = 3, inherit.aes = FALSE) +
    facet_wrap(vars(metric), ncol = 3, scales = "free_x") +
    labs(
      title    = "Gene-Segment Diversity Metrics vs True Clonotype Diversity",
      subtitle = "Red dotted line: 1:1 reference | Each facet = different gene-segment proxy",
      x        = "H (gene-segment diversity proxy)",
      y        = expression(H[clonotype]~"(true CDR3-level diversity)"),
      color    = group_var,
      caption  = paste("TRBD note: only ~17% of cells have D-gene assignment in J1568.",
                       "VDJ uses D-assigned cells only.")
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold", size = 10)
    )
}

# Internal helper: detect which columns in a wide table are the group_vars
# given that we only know one of them.
group_vars_from_wide <- function(wide, group_var) {
  # Return all columns that appear to be group keys (non-numeric, present in table)
  potential <- c("sample_id", "subset", "condition", "patient_id", group_var)
  intersect(potential, colnames(wide))
}

#' Bar chart of Spearman rho for each metric against true clonotype diversity
#'
#' Provides a quick visual ranking of which gene-segment combination best
#' tracks true clonotype diversity.
#'
#' @param correlations  The $correlations data frame from run_multisegment_analysis().
#' @param metric_order  Optional character vector giving display order of metrics.
#' @return  ggplot2 object.
plot_correlation_bars <- function(correlations,
                                   metric_order = c("TRBV", "TRBJ", "TRBD", "VJ", "VDJ")) {

  # Order metrics; only keep those present in data
  metric_order <- intersect(metric_order, correlations$metric)
  # Add any remaining metrics not in metric_order
  extra <- setdiff(correlations$metric, metric_order)
  metric_order <- c(metric_order, extra)

  df <- correlations %>%
    filter(metric %in% metric_order) %>%
    mutate(
      metric     = factor(metric, levels = rev(metric_order)),
      sig        = spearman_p < 0.05,
      rho_label  = paste0("rho = ", sprintf("%.2f", spearman_rho),
                          if_else(sig, "*", ""))
    )

  ggplot(df, aes(x = spearman_rho, y = metric, fill = sig)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = rho_label),
              hjust = if_else(df$spearman_rho >= 0, -0.05, 1.05),
              size = 3.5) +
    geom_vline(xintercept = 0, color = "grey40") +
    scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#92C5DE"),
                      labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05"),
                      name   = NULL) +
    scale_x_continuous(limits = c(min(c(df$spearman_rho, 0)) - 0.15,
                                   max(c(df$spearman_rho, 0)) + 0.25)) +
    labs(
      title    = "Spearman Correlation with True Clonotype Diversity (H_clonotype)",
      subtitle = "Higher rho = gene-segment proxy tracks true diversity better",
      x        = "Spearman rho",
      y        = "Gene-segment metric",
      caption  = "* p < 0.05 | TRBD: low power due to 83% missing D assignment in J1568"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")
}

#' Scatterplot matrix: all gene-segment H metrics vs each other and vs H_clonotype
#'
#' Useful for spotting collinearity between segment metrics.
#'
#' @param ms_results  Output of run_multisegment_analysis().
#' @return  ggplot2 object (pairs-style panel via patchwork).
plot_multisegment_pairs <- function(ms_results) {

  wide <- ms_results$metrics
  h_cols <- grep("^H_(TRBV|TRBJ|TRBD|VJ|VDJ|clonotype)$",
                 colnames(wide), value = TRUE)

  if (length(h_cols) < 2) {
    warning("Fewer than 2 H_ columns found; skipping pairs plot.")
    return(invisible(NULL))
  }

  # Build all unique pairs
  pairs <- combn(h_cols, 2, simplify = FALSE)

  plots <- lapply(pairs, function(p) {
    xc <- p[1]; yc <- p[2]
    xl <- sub("H_", "", xc); yl <- sub("H_", "", yc)
    sp <- tryCatch(
      cor.test(wide[[xc]], wide[[yc]], method = "spearman", exact = FALSE),
      error = function(e) NULL
    )
    sub_text <- if (!is.null(sp))
      paste0("rho = ", round(sp$estimate, 2))
    else ""

    ggplot(wide, aes(x = .data[[xc]], y = .data[[yc]])) +
      geom_point(size = 2, alpha = 0.8, color = "#2166AC") +
      geom_smooth(method = "lm", se = FALSE, color = "grey50",
                  linetype = "dashed", linewidth = 0.6) +
      labs(title = paste(xl, "vs", yl), subtitle = sub_text,
           x = paste("H", xl), y = paste("H", yl)) +
      theme_bw(base_size = 9)
  })

  wrap_plots(plots, ncol = ceiling(sqrt(length(plots))))
}

#' Print multi-segment correlation summary to console
#'
#' @param ms_results  Output of run_multisegment_analysis().
print_multisegment_summary <- function(ms_results) {
  cat("\n=== Multi-Segment Diversity: Correlation with True Clonotype Diversity ===\n")
  cat("(Spearman rho of H_[metric] vs H_clonotype)\n\n")
  print(ms_results$correlations %>%
          select(metric, n_groups, spearman_rho, spearman_p,
                 pearson_r, pearson_p) %>%
          arrange(desc(abs(spearman_rho))),
        row.names = FALSE)
  cat("\nKey:\n")
  cat("  TRBV  — V-gene usage only (NanoString-style, 44 genes)\n")
  cat("  TRBJ  — J-gene usage only (13 genes)\n")
  cat("  TRBD  — D-gene usage only (2 genes; 83% missing in J1568)\n")
  cat("  VJ    — V+J gene combination (440 unique combos in J1568)\n")
  cat("  VDJ   — V+D+J gene combination (687 unique; D-assigned cells only)\n")
  cat("  True  — CDR3 clonotype diversity (2866 unique; ground truth)\n\n")
  invisible(NULL)
}


# =============================================================================
# SECTION 11: GENE EXPRESSION PROGRAM SCORING
# =============================================================================
#
# PURPOSE:
#   Score each cell for five canonical T cell functional programs, then
#   average scores per sample × subset group and test whether those group-level
#   program scores correlate with true clonotype diversity — independently of
#   and in combination with TRBV/VJ diversity metrics.
#
# RATIONALE:
#   Clonal expansion leaves transcriptional footprints: expanded effector clones
#   upregulate cytotoxic, exhaustion, and proliferation programs. If per-sample
#   mean program scores track clonotype diversity, they carry complementary
#   information to V-gene usage. With n=8 samples this is exploratory only;
#   all results are reported as univariate Spearman rho (no multivariate model).
#
# PROGRAMS:
#   Proliferation — MKI67, TOP2A, PCNA, MCM2-6, STMN1, TYMS, CCNB1, CDK1, ...
#   Cytotoxic     — GZMB, GZMK, GZMH, PRF1, GNLY, NKG7, IFNG, CCL4, CCL5, ...
#   Exhaustion    — PDCD1, HAVCR2, TIGIT, LAG3, TOX, ENTPD1, CTLA4, CXCL13, ...
#   Activation    — CD69, ICOS, IL2RA, HLA-DRA, IRF4, JUNB, TNFRSF9, CD44, ...
#   Memory        — TCF7, CCR7, IL7R, SELL, LEF1, KLF2, S1PR1, FOXO1, BCL6, ...
#
# SCORING METHOD:
#   1. Library-size normalise: counts / colSums * 1e4 (CP10K)
#   2. Log1p transform
#   3. Per-cell program score = mean of log-normalised expression across
#      detected program genes (genes with ≥1 count anywhere in the matrix)
#   This is equivalent to a simplified Seurat AddModuleScore without background
#   correction. Use AddModuleScore for publication; this function is designed to
#   run without Seurat installed.
# =============================================================================

#' Return the canonical T cell gene program lists
#'
#' @return  Named list of character vectors; each element is a program with
#'          its constituent gene symbols.
define_tcell_programs <- function() {
  list(
    Proliferation = c("MKI67","TOP2A","PCNA","MCM2","MCM3","MCM4","MCM5","MCM6",
                      "STMN1","TYMS","CCNB1","CDK1","BIRC5","AURKB","UBE2C","PTTG1"),
    Cytotoxic     = c("GZMB","GZMK","GZMH","PRF1","GNLY","NKG7","IFNG",
                      "CCL4","CCL5","FGFBP2","CX3CR1","KLRG1"),
    Exhaustion    = c("PDCD1","HAVCR2","TIGIT","LAG3","TOX","ENTPD1",
                      "CTLA4","BATF","NR4A1","PRDM1","CXCL13","VCAM1"),
    Activation    = c("CD69","ICOS","IL2RA","HLA-DRA","IRF4","NFKBID",
                      "JUNB","DUSP6","TNFRSF9","TNFRSF4","CD44"),
    Memory        = c("TCF7","CCR7","IL7R","SELL","LEF1","KLF2",
                      "S1PR1","FOXO1","BCL6","CXCR5","ID3")
  )
}

#' Score T cell gene programs from a cellranger aggr expression matrix
#'
#' Reads the aggr sparse matrix, log-normalises, and computes a per-cell mean
#' score for each program in `programs`. Returns one row per cell.
#'
#' @param aggr_dir   Path to cellranger aggr output directory.
#' @param programs   Named list of gene symbol vectors (default: define_tcell_programs()).
#' @param barcodes   Optional character vector of barcodes to restrict scoring
#'                   to (e.g., only VDJ-matched T cells). NULL = score all cells.
#' @return  data.table with columns: barcode_aggr, <program_name> (score), ...
#'          plus <program_name>_n_genes (number of detected genes in program).
score_gene_programs <- function(aggr_dir,
                                 programs = define_tcell_programs(),
                                 barcodes = NULL) {

  mtx_file  <- file.path(aggr_dir, "matrix.mtx")
  feat_file <- file.path(aggr_dir, "features.tsv")
  bc_file   <- file.path(aggr_dir, "barcodes.tsv")

  message("  Reading expression matrix for program scoring...")
  mat  <- Matrix::readMM(mtx_file)   # features x cells (sparse)
  feat <- fread(feat_file, header = FALSE,
                col.names = c("ensembl","symbol","type"), showProgress = FALSE)
  feat[, row_idx := .I]
  bc_all <- fread(bc_file, header = FALSE, col.names = "barcode_aggr",
                  showProgress = FALSE)

  # Identify which cells to score
  if (!is.null(barcodes)) {
    cell_mask <- bc_all$barcode_aggr %in% barcodes
  } else {
    cell_mask <- rep(TRUE, nrow(bc_all))
  }
  mat_sub <- mat[, cell_mask, drop = FALSE]

  # Library-size normalisation → log1p
  # CP10K: counts per 10,000
  lib_sizes <- Matrix::colSums(mat_sub)
  lib_sizes[lib_sizes == 0] <- 1  # avoid division by zero
  mat_norm <- Matrix::t(Matrix::t(mat_sub) / lib_sizes) * 1e4
  mat_norm <- log1p(mat_norm)

  # Collect all unique program genes present in the matrix
  all_prog_genes <- unique(unlist(programs))
  gene_rows <- feat[symbol %in% all_prog_genes, .(symbol, row_idx)]

  # Extract dense sub-matrix for program genes only
  prog_dense <- as.matrix(mat_norm[gene_rows$row_idx, , drop = FALSE])
  rownames(prog_dense) <- gene_rows$symbol

  # Score each program per cell
  score_list <- lapply(names(programs), function(prog_name) {
    genes_in_prog   <- programs[[prog_name]]
    genes_available <- intersect(genes_in_prog, rownames(prog_dense))
    n_avail <- length(genes_available)

    if (n_avail == 0) {
      message("    WARNING: no genes found for program '", prog_name, "'")
      return(list(score = rep(NA_real_, sum(cell_mask)), n = 0L))
    }

    if (n_avail < length(genes_in_prog)) {
      missing <- setdiff(genes_in_prog, genes_available)
      message("    '", prog_name, "': using ", n_avail, "/", length(genes_in_prog),
              " genes (missing: ", paste(missing, collapse = ", "), ")")
    }

    sub_mat <- prog_dense[genes_available, , drop = FALSE]
    # Per-cell mean: only count genes with > 0 expression in that cell
    # (avoids diluting score with structural zeros)
    score <- colMeans(sub_mat)          # simple mean across all program genes
    list(score = score, n = n_avail)
  })
  names(score_list) <- names(programs)

  # Assemble result data.table
  result <- data.table(
    barcode_aggr = bc_all$barcode_aggr[cell_mask]
  )
  for (prog_name in names(score_list)) {
    result[, (prog_name)                          := score_list[[prog_name]]$score]
    result[, paste0(prog_name, "_n_genes") := score_list[[prog_name]]$n]
  }

  message("  Scored ", nrow(result), " cells across ",
          length(programs), " programs")
  result[]
}

#' Compute mean program scores per group (sample × subset)
#'
#' Averages per-cell program scores within each group defined by group_vars.
#' Returns one row per group with mean and fraction-positive for each program.
#'
#' @param cell_scores   data.table from score_gene_programs().
#' @param cell_meta     data.frame with barcode_aggr and group_vars columns.
#' @param group_vars    Character vector of grouping columns.
#' @param programs      Named list of gene programs (for column names).
#' @param min_cells     Minimum cells per group (default 20).
#' @return  data.table with one row per group and aggregated program scores.
aggregate_program_scores <- function(cell_scores,
                                      cell_meta,
                                      group_vars   = "sample_id",
                                      programs     = define_tcell_programs(),
                                      min_cells    = 20) {

  prog_names <- names(programs)

  # Join scores with metadata
  meta_dt <- as.data.table(cell_meta)[, c("barcode_aggr", group_vars), with = FALSE]
  joined  <- merge(meta_dt, cell_scores[, c("barcode_aggr", prog_names), with = FALSE],
                   by = "barcode_aggr", all.x = FALSE)  # inner join — drop unscored cells

  # Aggregate: mean score + fraction of cells with score > 0 per group
  agg <- joined[, {
    n_cells <- .N
    result  <- list(n_cells_scored = n_cells)
    for (p in prog_names) {
      vals <- .SD[[p]]
      result[[paste0("mean_", p)]]   <- mean(vals, na.rm = TRUE)
      result[[paste0("frac_pos_", p)]] <- mean(vals > 0, na.rm = TRUE)
    }
    result
  }, by = group_vars, .SDcols = prog_names]

  agg[n_cells_scored >= min_cells]
}

#' Correlate group-level program scores with diversity metrics
#'
#' For each program, computes Spearman rho with each diversity metric column
#' (H_clonotype, H_TRBV, H_VJ, etc.) across groups.
#'
#' @param program_agg   data.table from aggregate_program_scores().
#' @param diversity_tbl data.frame with diversity metrics (joined by group_vars).
#' @param group_vars    Character vector of join keys.
#' @param diversity_cols Character vector of diversity columns to correlate with.
#' @param programs       Named list of programs (for naming).
#' @return  data.frame with one row per (program × diversity_metric) pair:
#'          program, diversity_metric, spearman_rho, spearman_p, pearson_r,
#'          pearson_p, n_groups.
correlate_programs_with_diversity <- function(program_agg,
                                               diversity_tbl,
                                               group_vars     = "sample_id",
                                               diversity_cols = c("H_clonotype",
                                                                   "H_TRBV",
                                                                   "H_VJ"),
                                               programs       = define_tcell_programs()) {

  prog_score_cols <- paste0("mean_", names(programs))

  # Join
  joined <- merge(
    as.data.frame(program_agg),
    as.data.frame(diversity_tbl)[, c(group_vars, diversity_cols), drop = FALSE],
    by = group_vars
  )

  if (nrow(joined) < 3) {
    warning("Fewer than 3 groups after joining; correlations not computed.")
    return(NULL)
  }

  # Compute all pairwise correlations
  results <- purrr::map_dfr(prog_score_cols, function(pcol) {
    prog_label <- sub("^mean_", "", pcol)
    purrr::map_dfr(diversity_cols, function(dcol) {
      x <- joined[[pcol]]; y <- joined[[dcol]]
      ok <- !is.na(x) & !is.na(y)
      if (sum(ok) < 3) {
        return(data.frame(program = prog_label, diversity_metric = dcol,
                          spearman_rho = NA_real_, spearman_p = NA_real_,
                          pearson_r = NA_real_,    pearson_p  = NA_real_,
                          n_groups = sum(ok)))
      }
      sp <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
      pe <- cor.test(x[ok], y[ok], method = "pearson")
      data.frame(
        program          = prog_label,
        diversity_metric = dcol,
        spearman_rho     = round(sp$estimate, 4),
        spearman_p       = round(sp$p.value,  4),
        pearson_r        = round(pe$estimate, 4),
        pearson_p        = round(pe$p.value,  4),
        n_groups         = sum(ok)
      )
    })
  })

  results %>% arrange(diversity_metric, desc(abs(spearman_rho)))
}

#' Run the full program-diversity correlation pipeline
#'
#' Convenience wrapper: scores programs, aggregates per group, correlates
#' with diversity metrics, and returns everything needed for plotting.
#'
#' @param aggr_dir       cellranger aggr directory.
#' @param cell_meta      Data frame with barcode_aggr + group_vars columns.
#'                       Typically the VDJ metadata with subset annotations.
#' @param diversity_tbl  Joined diversity table (output of compare_diversity_metrics()$table
#'                       or multisegment metrics).
#' @param group_vars     Grouping columns (e.g. c("sample_id","subset")).
#' @param programs       Named list of gene programs.
#' @param min_cells      Minimum cells per group for aggregation.
#' @param diversity_cols Diversity metric columns to correlate with.
#' @return  List: $cell_scores, $group_scores, $correlations, $joined.
run_program_diversity_analysis <- function(aggr_dir,
                                            cell_meta,
                                            diversity_tbl,
                                            group_vars     = "sample_id",
                                            programs       = define_tcell_programs(),
                                            min_cells      = 20,
                                            diversity_cols = c("H_clonotype",
                                                               "H_TRBV", "H_VJ")) {

  message("=== Gene Program × Diversity Correlation Analysis ===")

  # Score programs per cell (restrict to cells in cell_meta for efficiency)
  message("\n[1/3] Scoring gene programs per cell...")
  cell_scores <- score_gene_programs(
    aggr_dir  = aggr_dir,
    programs  = programs,
    barcodes  = cell_meta$barcode_aggr
  )

  # Aggregate to group level
  message("\n[2/3] Aggregating scores per group...")
  group_scores <- aggregate_program_scores(
    cell_scores = cell_scores,
    cell_meta   = cell_meta,
    group_vars  = group_vars,
    programs    = programs,
    min_cells   = min_cells
  )
  message("  -> ", nrow(group_scores), " groups with sufficient cells")

  # Correlate with diversity metrics
  message("\n[3/3] Computing correlations with diversity metrics...")
  diversity_cols_available <- intersect(diversity_cols, colnames(diversity_tbl))
  corr <- correlate_programs_with_diversity(
    program_agg    = group_scores,
    diversity_tbl  = diversity_tbl,
    group_vars     = group_vars,
    diversity_cols = diversity_cols_available,
    programs       = programs
  )

  # Full joined table (program scores + diversity metrics)
  joined <- merge(
    as.data.frame(group_scores),
    as.data.frame(diversity_tbl),
    by = group_vars
  )

  message("\nDone.")
  list(
    cell_scores  = cell_scores,
    group_scores = group_scores,
    correlations = corr,
    joined       = joined
  )
}


# =============================================================================
# SECTION 12: PROGRAM × DIVERSITY PLOTTING
# =============================================================================

#' Heatmap of Spearman rho: gene programs (rows) × diversity metrics (columns)
#'
#' @param correlations  Output of correlate_programs_with_diversity().
#' @param metric_focus  Which diversity_metric values to include. NULL = all.
#' @param sig_threshold p-value threshold for asterisk annotation.
#' @return  Invisible NULL (pheatmap renders directly to device).
plot_program_correlation_heatmap <- function(correlations,
                                              metric_focus  = c("H_clonotype",
                                                                 "H_TRBV","H_VJ"),
                                              sig_threshold = 0.05) {

  df <- correlations
  if (!is.null(metric_focus)) {
    df <- df %>% filter(diversity_metric %in% metric_focus)
  }

  # Pretty metric labels
  metric_labels <- c(
    H_clonotype = "H clonotype\n(true CDR3 diversity)",
    H_TRBV      = "H TRBV\n(V-gene proxy)",
    H_VJ        = "H VJ\n(VJ combo proxy)"
  )

  df <- df %>%
    mutate(
      metric_label = if_else(diversity_metric %in% names(metric_labels),
                              metric_labels[diversity_metric],
                              diversity_metric),
      sig_label    = case_when(
        is.na(spearman_p)          ~ "",
        spearman_p < 0.001         ~ "***",
        spearman_p < 0.01          ~ "**",
        spearman_p < sig_threshold ~ "*",
        TRUE                       ~ ""
      ),
      cell_label = paste0(sprintf("%.2f", spearman_rho), sig_label)
    )

  # Build matrix
  rho_mat <- df %>%
    select(program, metric_label, spearman_rho) %>%
    pivot_wider(names_from = metric_label, values_from = spearman_rho) %>%
    tibble::column_to_rownames("program") %>%
    as.matrix()

  ann_mat <- df %>%
    select(program, metric_label, cell_label) %>%
    pivot_wider(names_from = metric_label, values_from = cell_label) %>%
    tibble::column_to_rownames("program") %>%
    as.matrix()

  # Reorder columns to match metric_focus order
  col_order <- intersect(metric_labels[metric_focus], colnames(rho_mat))
  rho_mat <- rho_mat[, col_order, drop = FALSE]
  ann_mat <- ann_mat[, col_order, drop = FALSE]

  pheatmap::pheatmap(
    rho_mat,
    display_numbers = ann_mat,
    number_color    = "black",
    fontsize_number = 11,
    color           = colorRampPalette(c("#D6604D","white","#2166AC"))(101),
    breaks          = seq(-1, 1, length.out = 102),
    cluster_rows    = TRUE,
    cluster_cols    = FALSE,
    border_color    = "grey90",
    main            = "Spearman rho: Gene Programs vs Diversity Metrics\n(* p<0.05  ** p<0.01  *** p<0.001)",
    fontsize_row    = 11,
    fontsize_col    = 10,
    angle_col       = 0
  )
  invisible(NULL)
}

#' Scatter panel: each gene program score vs H_clonotype (one facet per program)
#'
#' @param joined       $joined table from run_program_diversity_analysis().
#' @param y_col        Diversity metric for y-axis (default "H_clonotype").
#' @param color_var    Column name for point colour (e.g. "subset").
#' @param label_var    Column name for text labels (e.g. "sample_id").
#' @param programs     Named list of programs (for facet ordering).
#' @return  ggplot2 object.
plot_program_scatter_panel <- function(joined,
                                        y_col     = "H_clonotype",
                                        color_var = "subset",
                                        label_var = "sample_id",
                                        programs  = define_tcell_programs()) {

  prog_cols  <- paste0("mean_", names(programs))
  prog_cols  <- intersect(prog_cols, colnames(joined))
  prog_names <- sub("^mean_", "", prog_cols)

  # Correlation labels per facet
  rho_labels <- purrr::map_dfr(prog_cols, function(pc) {
    label <- sub("^mean_", "", pc)
    x <- joined[[pc]]; y <- joined[[y_col]]
    ok <- !is.na(x) & !is.na(y)
    if (sum(ok) < 3) return(data.frame(program=label, rho_label="n<3"))
    sp <- cor.test(x[ok], y[ok], method="spearman", exact=FALSE)
    data.frame(
      program   = label,
      rho_label = paste0("rho=", sprintf("%.2f", sp$estimate),
                          if_else(sp$p.value < 0.05, "*", ""),
                          "\np=", sprintf("%.3f", sp$p.value))
    )
  })

  # Long format
  plot_df <- purrr::map_dfr(prog_cols, function(pc) {
    label <- sub("^mean_", "", pc)
    joined %>%
      select(all_of(c(group_vars_from_wide(joined, color_var),
                      pc, y_col))) %>%
      rename(program_score = all_of(pc),
             diversity     = all_of(y_col)) %>%
      mutate(program = label)
  }) %>%
    mutate(program = factor(program, levels = names(programs)))

  # Colour palette
  n_colors <- if (color_var %in% colnames(plot_df)) n_distinct(plot_df[[color_var]]) else 1
  colors   <- if (n_colors == 2) c("#2166AC","#D6604D") else .trbv_palette(n_colors)

  p <- ggplot(plot_df, aes(x = program_score, y = diversity)) +
    geom_point(aes(color = if (color_var %in% colnames(plot_df))
                              .data[[color_var]] else "black"),
               size = 2.5, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, color = "grey40",
                linetype = "dashed", linewidth = 0.7) +
    geom_text_repel(aes(label = if (label_var %in% colnames(plot_df))
                                   .data[[label_var]] else ""),
                    size = 2.5, max.overlaps = 10) +
    geom_text(data = rho_labels,
              aes(label = rho_label),
              x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2,
              size = 3, inherit.aes = FALSE) +
    facet_wrap(vars(program), scales = "free_x", ncol = 3) +
    scale_color_manual(values = colors, name = color_var) +
    labs(
      title    = paste("Gene Program Scores vs", y_col),
      subtitle = "Mean log-normalised program expression per sample × subset group\n* p < 0.05 (Spearman)",
      x        = "Mean program score (log-normalised CP10K)",
      y        = y_col
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          strip.text = element_text(face = "bold"))
  p
}

#' Ranked bar chart: all programs and diversity proxy metrics, ranked by |rho| with H_clonotype
#'
#' Produces a single plot comparing programs and gene-segment metrics side by side.
#'
#' @param program_corr   Correlation table from correlate_programs_with_diversity().
#' @param segment_corr   Correlation table from run_multisegment_analysis()$correlations.
#' @param target_metric  Which diversity column to rank against (default "H_clonotype").
#' @return  ggplot2 object.
plot_combined_ranking <- function(program_corr,
                                   segment_corr,
                                   target_metric = "H_clonotype") {

  prog_rows <- program_corr %>%
    filter(diversity_metric == target_metric, !is.na(spearman_rho)) %>%
    transmute(feature = program, spearman_rho, spearman_p,
              feature_type = "Gene program")

  seg_rows <- segment_corr %>%
    filter(!is.na(spearman_rho)) %>%
    transmute(feature = metric, spearman_rho, spearman_p,
              feature_type = "Gene-segment metric")

  combined <- bind_rows(prog_rows, seg_rows) %>%
    mutate(
      sig   = !is.na(spearman_p) & spearman_p < 0.05,
      label = paste0(sprintf("%.2f", spearman_rho),
                     if_else(sig, "*", ""))
    ) %>%
    arrange(desc(spearman_rho)) %>%
    mutate(feature = factor(feature, levels = rev(feature)))

  ggplot(combined, aes(x = spearman_rho, y = feature,
                        fill = feature_type, alpha = sig)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = label),
              hjust = if_else(combined$spearman_rho >= 0, -0.08, 1.08),
              size = 3.2) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.5) +
    scale_fill_manual(values = c("Gene program"       = "#74ADD1",
                                  "Gene-segment metric" = "#D6604D"),
                       name = NULL) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.45),
                        guide = "none") +
    scale_x_continuous(
      limits = c(min(c(combined$spearman_rho, 0)) - 0.15,
                  max(c(combined$spearman_rho, 0)) + 0.25)
    ) +
    labs(
      title    = paste("All Features Ranked by Correlation with", target_metric),
      subtitle = "* p < 0.05 | Faded bars = non-significant",
      x        = "Spearman rho",
      y        = NULL,
      caption  = "Gene-segment metrics: TRBV, TRBJ, VJ (computed from TCR-seq).\nGene programs: mean log-normalised expression across program genes per group."
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")
}

#' Scatter plot of combined model: H_VJ + top program score vs H_clonotype
#'
#' Fits a simple rank-based additive model: rank(H_VJ) + rank(top_program).
#' With n=8 groups this is purely exploratory; shown as a combined rank score.
#'
#' @param joined        $joined table from run_program_diversity_analysis().
#' @param vj_col        Column name for VJ diversity (default "H_VJ").
#' @param top_program   Name of the program to combine (e.g. "Cytotoxic").
#' @param color_var     Column name for point colour.
#' @return  ggplot2 object.
plot_combined_rank_score <- function(joined,
                                      vj_col      = "H_VJ",
                                      top_program = "Cytotoxic",
                                      color_var   = "subset") {

  prog_col <- paste0("mean_", top_program)
  if (!prog_col %in% colnames(joined)) stop("Column not found: ", prog_col)
  if (!vj_col %in% colnames(joined))  stop("Column not found: ", vj_col)

  df <- joined %>%
    filter(!is.na(.data[[vj_col]]), !is.na(.data[[prog_col]]),
           !is.na(H_clonotype)) %>%
    mutate(
      # Additive rank model — rank both features, sum, rescale to 0-1
      rank_vj      = rank(.data[[vj_col]]),
      rank_prog    = rank(.data[[prog_col]]),
      combined_rank = (rank_vj + rank_prog) / (2 * n()),
      label         = if (color_var %in% colnames(.)) .data[[color_var]] else ""
    )

  sp_vj   <- cor.test(df[[vj_col]],       df$H_clonotype, method="spearman", exact=FALSE)
  sp_prog <- cor.test(df[[prog_col]],      df$H_clonotype, method="spearman", exact=FALSE)
  sp_comb <- cor.test(df$combined_rank,    df$H_clonotype, method="spearman", exact=FALSE)

  ann <- data.frame(
    model = c(paste("H_VJ only"),
              paste(top_program, "only"),
              paste("Combined rank\n(H_VJ +", top_program, ")")),
    rho   = c(sp_vj$estimate, sp_prog$estimate, sp_comb$estimate),
    p     = c(sp_vj$p.value,  sp_prog$p.value,  sp_comb$p.value)
  ) %>%
    mutate(label = paste0("rho=", sprintf("%.2f", rho),
                           if_else(p < 0.05, "*", ""),
                           " (p=", sprintf("%.3f", p), ")"))

  cat("Combined rank model summary:\n")
  print(ann[, c("model","rho","p","label")])

  n_colors <- if (color_var %in% colnames(df)) n_distinct(df[[color_var]]) else 1
  colors   <- if (n_colors == 2) c("#2166AC","#D6604D") else .trbv_palette(n_colors)

  ggplot(df, aes(x = combined_rank, y = H_clonotype,
                  color = if (color_var %in% colnames(df))
                            .data[[color_var]] else "grey40")) +
    geom_point(size = 3, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, color = "grey40",
                linetype = "dashed", linewidth = 0.7) +
    geom_text_repel(aes(label = if ("sample_id" %in% colnames(df))
                                   sample_id else ""),
                    size = 3, max.overlaps = 15) +
    scale_color_manual(values = colors, name = color_var) +
    annotate("text", x = 0.05, y = max(df$H_clonotype, na.rm=TRUE),
             label = paste0("Combined: ", ann$label[3], "\n",
                             "H_VJ alone: ", ann$label[1], "\n",
                             top_program, " alone: ", ann$label[2]),
             hjust = 0, vjust = 1, size = 3.2, color = "grey20") +
    labs(
      title    = paste0("Combined Rank Score (H_VJ + ", top_program,
                         ") vs True Clonotype Diversity"),
      subtitle = "Rank-based additive combination — exploratory only (n=8)",
      x        = paste0("Combined rank score [rank(H_VJ) + rank(", top_program, ")] / 2n"),
      y        = expression(H[clonotype]~"(true CDR3 diversity)")
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")
}


# =============================================================================
# SECTION 12: CANONICAL TCR GENE LISTS AND FULL-CHAIN DIVERSITY
# =============================================================================
#
# Extends the TRBV-only proxy to all canonical human TCR gene segments:
# TRBV, TRBD, TRBJ (beta chain) and TRAV, TRAJ (alpha chain), plus paired
# alpha-beta combination metrics.
# =============================================================================

# ── Canonical gene lists (IMGT human) ────────────────────────────────────────

TCR_GENES <- list(

  TRBV = c(
    "TRBV1","TRBV2","TRBV3-1","TRBV4-1","TRBV4-2","TRBV4-3",
    "TRBV5-1","TRBV5-4","TRBV5-5","TRBV5-6","TRBV5-8",
    "TRBV6-1","TRBV6-2","TRBV6-3","TRBV6-4","TRBV6-5","TRBV6-6","TRBV6-7","TRBV6-8","TRBV6-9",
    "TRBV7-2","TRBV7-3","TRBV7-4","TRBV7-6","TRBV7-7","TRBV7-8","TRBV7-9",
    "TRBV9","TRBV10-1","TRBV10-2","TRBV10-3",
    "TRBV11-1","TRBV11-2","TRBV11-3",
    "TRBV12-1","TRBV12-2","TRBV12-3","TRBV12-4",
    "TRBV13","TRBV14","TRBV15","TRBV16","TRBV17","TRBV18","TRBV19","TRBV20-1",
    "TRBV24-1","TRBV25-1","TRBV27","TRBV28","TRBV29-1","TRBV30"
  ),

  TRBD = c("TRBD1","TRBD2"),

  TRBJ = c(
    "TRBJ1-1","TRBJ1-2","TRBJ1-3","TRBJ1-4","TRBJ1-5","TRBJ1-6",
    "TRBJ2-1","TRBJ2-2","TRBJ2-3","TRBJ2-4","TRBJ2-5","TRBJ2-6","TRBJ2-7"
  ),

  TRAV = c(
    "TRAV1-1","TRAV1-2","TRAV2","TRAV3","TRAV4","TRAV5","TRAV6","TRAV7",
    "TRAV8-1","TRAV8-2","TRAV8-3","TRAV8-4","TRAV9-1","TRAV9-2","TRAV10",
    "TRAV12-1","TRAV12-2","TRAV12-3","TRAV13-1","TRAV13-2","TRAV14DV4",
    "TRAV16","TRAV17","TRAV19","TRAV20","TRAV21","TRAV22","TRAV23DV6",
    "TRAV24","TRAV25","TRAV26-1","TRAV26-2","TRAV27","TRAV29DV5","TRAV30",
    "TRAV34","TRAV35","TRAV36DV7","TRAV38-1","TRAV38-2DV8","TRAV39","TRAV40","TRAV41"
  ),

  TRAJ = c(
    "TRAJ1","TRAJ2","TRAJ3","TRAJ4","TRAJ5","TRAJ6","TRAJ7","TRAJ8","TRAJ9","TRAJ10",
    "TRAJ11","TRAJ12","TRAJ13","TRAJ14","TRAJ15","TRAJ16","TRAJ17","TRAJ18","TRAJ19","TRAJ20",
    "TRAJ21","TRAJ22","TRAJ23","TRAJ24","TRAJ25","TRAJ26","TRAJ27","TRAJ28","TRAJ29","TRAJ30",
    "TRAJ31","TRAJ32","TRAJ33","TRAJ34","TRAJ35","TRAJ36","TRAJ37","TRAJ38","TRAJ39","TRAJ40",
    "TRAJ41","TRAJ42","TRAJ43","TRAJ44","TRAJ45","TRAJ46","TRAJ47","TRAJ48","TRAJ49","TRAJ50",
    "TRAJ51","TRAJ52","TRAJ53","TRAJ54","TRAJ55","TRAJ56","TRAJ57","TRAJ58","TRAJ59","TRAJ60",
    "TRAJ61"
  )
)

#' Standardise gene calls against canonical IMGT gene lists
#'
#' Strips allele suffixes (*01), trims whitespace, and optionally replaces
#' unrecognised calls with NA.
#'
#' @param x              Character vector of raw gene calls.
#' @param chain          One of "TRBV","TRBD","TRBJ","TRAV","TRAJ" (or NULL).
#' @param na_if_unknown  If TRUE, replace unrecognised calls with NA.
#' @return  Cleaned character vector.
standardise_gene_calls <- function(x, chain = NULL, na_if_unknown = FALSE) {
  x <- str_trim(as.character(x))
  x <- sub("\\*\\d+$", "", x)
  x[x %in% c("", "None", "none", "nan", "NaN", "NA", "N/A", "na")] <- NA_character_
  if (!is.null(chain) && na_if_unknown && chain %in% names(TCR_GENES)) {
    x[!is.na(x) & !x %in% TCR_GENES[[chain]]] <- NA_character_
  }
  x
}

#' Add standardised TRAV/TRAJ columns and all paired combo columns
#'
#' Adds: trav_gene_std, traj_gene_std, vj_alpha (TRAV:TRAJ),
#' vavb_combo (TRAV:TRBV), full_vj_combo (TRAV:TRAJ:TRBV:TRBJ).
#'
#' @param metadata       Data frame. Needs trav_gene and traj_gene columns.
#' @param na_if_unknown  If TRUE, replace off-list gene calls with NA.
#' @return  metadata with additional columns.
add_alpha_chain_combos <- function(metadata, na_if_unknown = FALSE) {
  for (m in c("trav_gene", "traj_gene")) {
    if (!m %in% colnames(metadata)) {
      warning("add_alpha_chain_combos: column '", m, "' not found; setting to NA.")
      metadata[[m]] <- NA_character_
    }
  }
  metadata <- metadata %>%
    mutate(
      trav_gene_std = standardise_gene_calls(trav_gene, "TRAV", na_if_unknown),
      traj_gene_std = standardise_gene_calls(traj_gene, "TRAJ", na_if_unknown),
      vj_alpha = if_else(
        !is.na(trav_gene_std) & !is.na(traj_gene_std),
        paste(trav_gene_std, traj_gene_std, sep = ":"),
        NA_character_
      ),
      vavb_combo = if_else(
        !is.na(trav_gene_std) & "trbv_gene" %in% colnames(.) &
          !is.na(trbv_gene) & trbv_gene != "",
        paste(trav_gene_std, trbv_gene, sep = ":"),
        NA_character_
      ),
      full_vj_combo = if_else(
        !is.na(trav_gene_std) & !is.na(traj_gene_std) &
          "trbv_gene" %in% colnames(.) & !is.na(trbv_gene) & trbv_gene != "" &
          "trbj_gene" %in% colnames(.) & !is.na(trbj_gene) & trbj_gene != "",
        paste(trav_gene_std, traj_gene_std, trbv_gene, trbj_gene, sep = ":"),
        NA_character_
      )
    )
  metadata
}

#' Run complete multi-segment TCR diversity analysis (beta + alpha + paired)
#'
#' Covers TRBV, TRBD, TRBJ, TRAV, TRAJ, VJ_beta, VJ_alpha, VA:VB, and the
#' full four-segment paired proxy. Correlates each with true clonotype diversity.
#'
#' @param metadata    Data frame after add_vdj_combos() and add_alpha_chain_combos().
#' @param group_vars  Grouping column(s).
#' @param min_cells   Minimum cells per group.
#' @param downsample_n Optional downsampling target.
#' @return  List with $metrics, $correlations, $clonotype.
run_full_tcr_analysis <- function(metadata,
                                   group_vars   = "sample_id",
                                   min_cells    = 30,
                                   downsample_n = 0) {

  if (!"vj_combo" %in% colnames(metadata))
    metadata <- add_vdj_combos(metadata)
  if (!"vj_alpha" %in% colnames(metadata))
    metadata <- add_alpha_chain_combos(metadata)

  # Standardise beta-chain calls
  for (col in c("trbv_gene","trbj_gene","trbd_gene")) {
    chain <- toupper(sub("_gene$","",col))
    if (col %in% colnames(metadata))
      metadata[[col]] <- standardise_gene_calls(metadata[[col]], chain)
  }

  message("=== Full TCR Multi-Segment Analysis (Beta + Alpha + Paired) ===")

  seg_info <- list(
    list(col="trbv_gene",     label="TRBV"),
    list(col="trbj_gene",     label="TRBJ"),
    list(col="trbd_gene",     label="TRBD"),
    list(col="trav_gene_std", label="TRAV"),
    list(col="traj_gene_std", label="TRAJ"),
    list(col="vj_combo",      label="VJ_beta"),
    list(col="vj_alpha",      label="VJ_alpha"),
    list(col="vavb_combo",    label="VA_VB"),
    list(col="full_vj_combo", label="Full_paired")
  )

  for (s in seg_info) {
    n <- if (s$col %in% colnames(metadata))
      sum(!is.na(metadata[[s$col]]) & metadata[[s$col]] != "") else 0L
    message(sprintf("  %-14s: %d / %d cells (%.1f%%)",
                    s$label, n, nrow(metadata), 100*n/nrow(metadata)))
  }

  clono      <- compute_clonotype_diversity(metadata, group_vars, min_cells, downsample_n)
  seg_tables <- lapply(seg_info, function(s) {
    if (!s$col %in% colnames(metadata)) return(NULL)
    n_valid <- sum(!is.na(metadata[[s$col]]) & metadata[[s$col]] != "")
    if (n_valid < min_cells) return(NULL)
    message("  Computing ", s$label, "...")
    compute_gene_segment_diversity(metadata, s$col, group_vars,
                                    min_cells, downsample_n, metric_label=s$label)
  })
  seg_tables <- Filter(Negate(is.null), seg_tables)

  metrics_wide <- Reduce(function(a,b) left_join(a, b, by=group_vars), seg_tables)
  metrics_wide <- left_join(metrics_wide, clono, by=group_vars)

  h_cols <- intersect(paste0("H_", vapply(seg_info,`[[`,character(1),"label")),
                      colnames(metrics_wide))

  correlations <- purrr::map_dfr(h_cols, function(hcol) {
    label <- sub("^H_","",hcol)
    x <- metrics_wide[[hcol]]; y <- metrics_wide$H_clonotype
    ok <- !is.na(x) & !is.na(y)
    if (sum(ok) < 3)
      return(data.frame(metric=label, n_groups=sum(ok),
                        spearman_rho=NA_real_, spearman_p=NA_real_,
                        pearson_r=NA_real_, pearson_p=NA_real_))
    sp <- cor.test(x[ok], y[ok], method="spearman", exact=FALSE)
    pe <- cor.test(x[ok], y[ok], method="pearson")
    data.frame(metric=label, n_groups=sum(ok),
               spearman_rho=round(sp$estimate,4), spearman_p=round(sp$p.value,4),
               pearson_r=round(pe$estimate,4),    pearson_p=round(pe$p.value,4))
  }) %>% arrange(desc(abs(spearman_rho)))

  message("\nCorrelation ranking (Spearman ρ with H_clonotype):")
  print(correlations %>% select(metric, n_groups, spearman_rho, spearman_p))

  list(metrics=metrics_wide, correlations=correlations, clonotype=clono)
}

#' Stacked barplot for any TCR gene segment column
#'
#' Works for TRAV, TRAJ, TRBV, TRBJ, TRBD, or any discrete gene call column.
#'
#' @param metadata    Data frame with gene_col and group_var columns.
#' @param gene_col    Column name of the gene calls.
#' @param group_var   Column for x-axis grouping.
#' @param top_n_genes Top N genes to show (remainder → "Other").
#' @param title       Optional plot title.
#' @return  ggplot2 object.
plot_gene_usage_bar <- function(metadata,
                                 gene_col,
                                 group_var   = "sample_id",
                                 top_n_genes = 20,
                                 title       = NULL) {
  top_genes <- metadata %>%
    filter(!is.na(.data[[gene_col]]), .data[[gene_col]] != "") %>%
    count(.data[[gene_col]], sort=TRUE) %>%
    slice_head(n=top_n_genes) %>%
    pull(1)

  plot_df <- metadata %>%
    filter(!is.na(.data[[gene_col]]), .data[[gene_col]] != "") %>%
    mutate(gene_plot = if_else(.data[[gene_col]] %in% top_genes,
                                .data[[gene_col]], "Other")) %>%
    count(.data[[group_var]], gene_plot) %>%
    group_by(.data[[group_var]]) %>%
    mutate(freq = n / sum(n)) %>%
    ungroup() %>%
    mutate(gene_plot = factor(gene_plot, levels=c(sort(top_genes),"Other")))

  n_col  <- length(unique(plot_df$gene_plot))
  colors <- c(.trbv_palette(n_col - 1), "grey70")
  names(colors) <- levels(plot_df$gene_plot)

  ggplot(plot_df, aes(x=.data[[group_var]], y=freq, fill=gene_plot)) +
    geom_bar(stat="identity", width=0.8, color="white", linewidth=0.2) +
    scale_fill_manual(values=colors, name=gene_col) +
    scale_y_continuous(labels=percent_format()) +
    labs(title    = title %||% paste(gene_col, "Usage per Sample"),
         subtitle = paste0("Top ", top_n_genes, " genes shown"),
         x = group_var, y = "Proportion of T cells") +
    theme_bw(base_size=12) +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          legend.key.size=unit(0.4,"cm"),
          legend.text=element_text(size=8))
}

#' Heatmap of gene usage proportions for any TCR segment
#'
#' @param metadata    Data frame with gene_col and group_var columns.
#' @param gene_col    Column name of the gene calls.
#' @param group_var   Column for heatmap columns.
#' @param top_n_genes Number of top genes to display.
#' @return  Invisible NULL (pheatmap renders directly).
plot_gene_usage_heatmap <- function(metadata,
                                     gene_col,
                                     group_var   = "sample_id",
                                     top_n_genes = 30) {
  top_genes <- metadata %>%
    filter(!is.na(.data[[gene_col]]), .data[[gene_col]] != "") %>%
    count(.data[[gene_col]], sort=TRUE) %>%
    slice_head(n=top_n_genes) %>%
    pull(1)

  heat_mat <- metadata %>%
    filter(.data[[gene_col]] %in% top_genes) %>%
    count(.data[[group_var]], .data[[gene_col]]) %>%
    group_by(.data[[group_var]]) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    select(-n) %>%
    pivot_wider(names_from=all_of(group_var),
                values_from=prop, values_fill=0) %>%
    tibble::column_to_rownames(gene_col) %>%
    as.matrix()

  pheatmap::pheatmap(
    heat_mat,
    color        = colorRampPalette(c("white","#2166AC","#053061"))(100),
    main         = paste(gene_col, "Usage Heatmap"),
    fontsize_row = 7, fontsize_col = 9,
    cluster_rows = TRUE, cluster_cols = TRUE, border_color = NA
  )
  invisible(NULL)
}
}


# =============================================================================
# SECTION 13: CD4/CD8 STRATIFICATION AND PSEUDOCLONOTYPE ANALYSIS
# =============================================================================
#
# OVERVIEW:
#   Part A — Subset detection and stratification
#     Detects cell type annotation columns in h5ad-derived metadata, maps them
#     to CD4 / CD8 / Treg / DN / Unknown, and enables all diversity analyses
#     to be run per sample × subset rather than per sample only.
#
#   Part B — Pseudoclonotype assignment and recapitulation metrics
#     For each gene-combination level (from coarse TRBV-only to the full
#     4-segment TRAV:TRAJ:TRBV:TRBJ paired proxy), cells sharing the same
#     gene-combo are assigned the same "pseudoclonotype" label. The analysis
#     then quantifies how well each level recapitulates true CDR3 clonotypes
#     using four complementary metrics:
#
#       Purity        — of cells sharing a pseudoclonotype, what fraction share
#                       the dominant real clonotype? (high = each pseudo maps to
#                       one real clone)
#       Collision rate — fraction of pseudoclonotypes containing >1 distinct real
#                       clonotype (high = proxy lumps unrelated clones together)
#       Split rate    — fraction of real clonotypes spanning >1 pseudoclonotype
#                       (high = proxy artificially fragments a real clone)
#       H-correlation — Spearman ρ of Shannon diversity(pseudoclonotype) vs
#                       Shannon diversity(real clonotype) across samples
#
#   All analyses are stratified by CD4/CD8 subset when annotation is present.
# =============================================================================


# ── PART A: CD4/CD8 SUBSET DETECTION AND STANDARDISATION ─────────────────────

# Common column name patterns that carry T cell subset annotations
.SUBSET_COL_CANDIDATES <- c(
  "cell_type","celltype","cell_type_l1","cell_type_l2","cell_type_l3",
  "subset","T_cell_type","T_cell_subset","lineage","CD4_CD8","cd4_cd8",
  "cluster_label","leiden_label","annotation","predicted.id","SingleR.labels",
  "majority_voting","sctype","azimuth","broad_cell_type","fine_cell_type"
)

# Regex patterns mapping annotation strings to canonical subsets
.SUBSET_MAP <- list(
  CD8  = c("CD8","cd8","cytotoxic","CTL","cytotox","TE","Tem","Tex",
           "exhausted","effector memory","effector CD8","CD8T"),
  CD4  = c("CD4","cd4","T helper","Th1","Th2","Th17","Tfh","T_helper",
           "helper","CD4T","naive CD4"),
  Treg = c("Treg","T_reg","regulatory","FOXP3","T regulatory"),
  DN   = c("double negative","DN_T","double-negative","DN T","gamma.delta","gd","gdT")
)

#' Auto-detect the best T cell subset annotation column in a metadata data frame
#'
#' Checks candidate column names, then for each candidate verifies that at
#' least 5% of non-NA values match a CD4 or CD8 pattern. Returns the column
#' with the highest fraction of mappable values.
#'
#' @param metadata  Data frame (from load_h5ad_metadata).
#' @return  Column name (character) or NULL if none found.
detect_subset_column <- function(metadata) {
  cols <- colnames(metadata)

  # Exact or substring match against candidates
  hits <- intersect(.SUBSET_COL_CANDIDATES, cols)
  if (length(hits) == 0) {
    # Fall back to fuzzy match
    hits <- cols[grep(paste(.SUBSET_COL_CANDIDATES, collapse="|"),
                      cols, ignore.case=TRUE)]
  }
  if (length(hits) == 0) return(NULL)

  cd4_pat <- paste(.SUBSET_MAP$CD4,  collapse="|")
  cd8_pat <- paste(.SUBSET_MAP$CD8,  collapse="|")

  scores <- vapply(hits, function(col) {
    vals <- as.character(metadata[[col]])
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) return(0)
    mean(grepl(cd4_pat, vals, ignore.case=TRUE) |
         grepl(cd8_pat, vals, ignore.case=TRUE))
  }, numeric(1))

  best <- hits[which.max(scores)]
  if (max(scores) < 0.05) return(NULL)  # < 5% mappable → not useful
  best
}

#' Map raw annotation strings to canonical T cell subsets
#'
#' Returns one of: "CD8", "CD4", "Treg", "DN", "Non_T", "Unknown".
#'
#' @param x   Character vector of raw cell type annotations.
#' @return  Character vector of canonical labels.
standardise_subset <- function(x) {
  x <- as.character(x)
  out <- rep("Unknown", length(x))
  out[is.na(x) | x %in% c("","None","nan","NaN")] <- "Unknown"

  for (label in c("Treg","CD8","CD4","DN")) {   # Treg before CD4 (subset)
    pat <- paste(.SUBSET_MAP[[label]], collapse="|")
    out[grepl(pat, x, ignore.case=TRUE)] <- label
  }
  out[grepl("non.?t|nonT|non_T|B cell|NK|Mono|DC|macro|plasma|mast|neutro",
            x, ignore.case=TRUE)] <- "Non_T"
  out
}

#' Add a standardised 'subset' column to TCR metadata
#'
#' Detects the best annotation column and maps values to CD4/CD8/Treg/DN.
#' Optionally removes Non_T and Unknown cells.
#'
#' @param metadata       Data frame with cell-level metadata.
#' @param subset_col     Override: exact column name to use (skips auto-detect).
#' @param keep_nont      If FALSE (default), remove Non_T cells.
#' @param keep_unknown   If FALSE (default), remove Unknown cells.
#' @return  metadata with 'subset' column added. Reports mapping counts.
add_subset_from_annotation <- function(metadata,
                                        subset_col   = NULL,
                                        keep_nont    = FALSE,
                                        keep_unknown = FALSE) {

  if (is.null(subset_col)) subset_col <- detect_subset_column(metadata)

  if (is.null(subset_col)) {
    warning("No subset annotation column detected. ",
            "Set subset_col= manually or use annotate_tcell_subsets_from_aggr().\n",
            "Returning metadata with subset = 'Unknown' for all cells.")
    metadata$subset <- "Unknown"
    return(metadata)
  }

  message("  Using '", subset_col, "' for T cell subset annotation.")
  metadata$subset <- standardise_subset(metadata[[subset_col]])

  tab <- table(metadata$subset)
  message("  Subset mapping:")
  for (nm in names(tab))
    message(sprintf("    %-10s: %d (%.1f%%)", nm, tab[nm], 100*tab[nm]/nrow(metadata)))

  if (!keep_nont)    metadata <- metadata[metadata$subset != "Non_T",    ]
  if (!keep_unknown) metadata <- metadata[metadata$subset != "Unknown",  ]

  metadata
}


# ── PART B: PSEUDOCLONOTYPE ASSIGNMENT AND RECAPITULATION METRICS ─────────────

# Gene-combo levels in increasing order of resolution
PSEUDO_LEVELS <- list(
  list(col="trbv_gene",     label="TRBV"),
  list(col="trbj_gene",     label="TRBJ"),
  list(col="vj_combo",      label="VJ_beta"),
  list(col="trav_gene_std", label="TRAV"),
  list(col="vj_alpha",      label="VJ_alpha"),
  list(col="vavb_combo",    label="VA_VB"),
  list(col="full_vj_combo", label="Full_paired")
)

#' Compute pseudoclonotype recapitulation metrics for one gene-combo level
#'
#' For each group (sample × subset), assigns pseudoclonotype labels based on
#' shared gene-combo strings, then computes purity, collision rate, split rate,
#' and diversity correlation vs real CDR3 clonotypes.
#'
#' @param metadata   Data frame. Must contain pseudo_col, clono_col,
#'                   and group_vars columns.
#' @param pseudo_col Column name of the gene-combo to use as pseudoclonotype
#'                   labels (e.g., "full_vj_combo", "vj_combo").
#' @param clono_col  Column name of real CDR3 clonotype labels. Default "clonotype_id".
#' @param group_vars Grouping column names (e.g., c("sample_id","subset")).
#' @param min_cells  Minimum cells per group (default 30).
#' @return  Data frame with one row per group and columns:
#'   n_cells, coverage, n_pseudo, n_real, pseudo_real_ratio,
#'   H_pseudo, H_real, H_norm_pseudo, H_norm_real,
#'   purity, collision_rate, split_rate, spearman_rho_H
compute_pseudoclonotype_metrics <- function(metadata,
                                             pseudo_col,
                                             clono_col  = "clonotype_id",
                                             group_vars = "sample_id",
                                             min_cells  = 30) {

  if (!pseudo_col %in% colnames(metadata))
    stop("Column '", pseudo_col, "' not found in metadata.")
  if (!clono_col  %in% colnames(metadata))
    stop("Column '", clono_col,  "' not found in metadata.")

  # Filter to cells with both assignments
  df <- metadata %>%
    filter(!is.na(.data[[pseudo_col]]), .data[[pseudo_col]] != "",
           !is.na(.data[[clono_col]]),  .data[[clono_col]]  != "")

  if (nrow(df) == 0) {
    warning("No cells with both '", pseudo_col, "' and '", clono_col, "' assigned.")
    return(tibble())
  }

  groups <- df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(n_cells_both = n(), .groups="drop") %>%
    filter(n_cells_both >= min_cells)

  if (nrow(groups) == 0) {
    warning("No groups pass min_cells threshold (", min_cells, ").")
    return(tibble())
  }

  purrr::map_dfr(seq_len(nrow(groups)), function(i) {
    # Filter to this group
    grp <- groups[i, , drop=FALSE]
    grp_df <- df
    for (gv in group_vars)
      grp_df <- grp_df[grp_df[[gv]] == grp[[gv]], ]

    n_cells    <- nrow(grp_df)
    n_cells_all <- sum(metadata[[group_vars[1]]] == grp[[group_vars[1]]])

    pseudo_vec <- grp_df[[pseudo_col]]
    real_vec   <- grp_df[[clono_col]]

    # Coverage = fraction of all cells in group with valid assignment
    coverage <- n_cells / max(n_cells_all, 1)

    # Shannon diversity
    pseudo_cts <- table(pseudo_vec)
    real_cts   <- table(real_vec)
    H_pseudo   <- shannon_diversity(as.numeric(pseudo_cts))
    H_real     <- shannon_diversity(as.numeric(real_cts))
    Hn_pseudo  <- normalized_shannon(as.numeric(pseudo_cts))
    Hn_real    <- normalized_shannon(as.numeric(real_cts))

    n_pseudo <- length(pseudo_cts)
    n_real   <- length(real_cts)

    # ── Purity ──────────────────────────────────────────────────────────────
    # For each pseudoclonotype, what fraction of its cells share the dominant
    # real clonotype?
    purity_per_pseudo <- grp_df %>%
      group_by(.data[[pseudo_col]]) %>%
      summarise(
        n_in_pseudo = n(),
        dom_frac = max(table(.data[[clono_col]])) / n(),
        .groups = "drop"
      )
    # Weighted purity (weight by pseudoclonotype size)
    purity <- with(purity_per_pseudo,
                   sum(n_in_pseudo * dom_frac) / sum(n_in_pseudo))

    # ── Collision rate ────────────────────────────────────────────────────────
    # Fraction of pseudoclonotypes that contain >1 distinct real clonotype
    n_collisions <- grp_df %>%
      group_by(.data[[pseudo_col]]) %>%
      summarise(n_real_in_pseudo = n_distinct(.data[[clono_col]]), .groups="drop") %>%
      summarise(n_coll = sum(n_real_in_pseudo > 1), n_total = n()) %>%
      { .$n_coll / .$n_total }
    collision_rate <- as.numeric(n_collisions)

    # ── Split rate ────────────────────────────────────────────────────────────
    # Fraction of real clonotypes that span >1 pseudoclonotype
    n_splits <- grp_df %>%
      group_by(.data[[clono_col]]) %>%
      summarise(n_pseudo_for_real = n_distinct(.data[[pseudo_col]]), .groups="drop") %>%
      summarise(n_split = sum(n_pseudo_for_real > 1), n_total = n()) %>%
      { .$n_split / .$n_total }
    split_rate <- as.numeric(n_splits)

    # Combine group keys + metrics
    bind_cols(
      grp,
      tibble(
        n_cells          = n_cells,
        coverage         = round(coverage, 3),
        n_pseudo         = n_pseudo,
        n_real           = n_real,
        pseudo_real_ratio= round(n_pseudo / max(n_real, 1), 3),
        H_pseudo         = round(H_pseudo,  4),
        H_real           = round(H_real,    4),
        H_norm_pseudo    = round(Hn_pseudo, 4),
        H_norm_real      = round(Hn_real,   4),
        purity           = round(purity,         3),
        collision_rate   = round(collision_rate,  3),
        split_rate       = round(split_rate,      3)
      )
    )
  })
}

#' Run pseudoclonotype recapitulation analysis across all gene-combo levels
#'
#' Calls compute_pseudoclonotype_metrics() for every level defined in
#' PSEUDO_LEVELS and computes a per-level correlation summary.
#'
#' @param metadata   Data frame with all combo columns added
#'                   (run add_vdj_combos() and add_alpha_chain_combos() first).
#' @param group_vars Grouping columns (e.g., c("sample_id","subset")).
#' @param min_cells  Minimum cells per group.
#' @return  List with:
#'   $per_level  — metrics data frame for each pseudo level (long format)
#'   $summary    — one row per level: mean purity, collision, split, n_groups
#'   $h_corr     — Spearman ρ of H_pseudo vs H_real per level
run_pseudoclonotype_analysis <- function(metadata,
                                          group_vars = "sample_id",
                                          min_cells  = 30) {

  message("=== Pseudoclonotype Recapitulation Analysis ===")

  results_list <- lapply(PSEUDO_LEVELS, function(lvl) {
    if (!lvl$col %in% colnames(metadata)) {
      message("  Skipping ", lvl$label, " (column not found)")
      return(NULL)
    }
    n_valid <- sum(!is.na(metadata[[lvl$col]]) & metadata[[lvl$col]] != "")
    if (n_valid < min_cells) {
      message("  Skipping ", lvl$label, " (< ", min_cells, " valid cells)")
      return(NULL)
    }
    message("  Computing metrics for ", lvl$label, " ...")
    tbl <- compute_pseudoclonotype_metrics(
      metadata, lvl$col, "clonotype_id", group_vars, min_cells
    )
    if (nrow(tbl) == 0) return(NULL)
    tbl$pseudo_level <- lvl$label
    tbl
  })
  results_list <- Filter(Negate(is.null), results_list)

  if (length(results_list) == 0) {
    warning("No pseudoclonotype results computed.")
    return(list(per_level=tibble(), summary=tibble(), h_corr=tibble()))
  }

  per_level <- bind_rows(results_list)

  # Summary per level
  summary_tbl <- per_level %>%
    group_by(pseudo_level) %>%
    summarise(
      n_groups       = n(),
      mean_purity    = round(mean(purity,         na.rm=TRUE), 3),
      sd_purity      = round(sd(purity,           na.rm=TRUE), 3),
      mean_collision = round(mean(collision_rate,  na.rm=TRUE), 3),
      mean_split     = round(mean(split_rate,      na.rm=TRUE), 3),
      mean_coverage  = round(mean(coverage,        na.rm=TRUE), 3),
      mean_n_pseudo  = round(mean(n_pseudo,        na.rm=TRUE), 1),
      mean_n_real    = round(mean(n_real,          na.rm=TRUE), 1),
      .groups = "drop"
    ) %>%
    # Preserve the resolution order
    mutate(pseudo_level = factor(pseudo_level,
                                  levels=vapply(PSEUDO_LEVELS,`[[`,character(1),"label"))) %>%
    arrange(pseudo_level)

  # H_pseudo vs H_real Spearman correlation per level
  h_corr <- per_level %>%
    group_by(pseudo_level) %>%
    summarise(
      n_groups     = n(),
      spearman_rho = {
        ok <- !is.na(H_pseudo) & !is.na(H_real) & n() >= 3
        if (sum(ok) < 3) NA_real_
        else round(cor(H_pseudo[ok], H_real[ok], method="spearman"), 4)
      },
      .groups = "drop"
    )

  message("\n--- Pseudoclonotype summary ---")
  print(summary_tbl %>% select(pseudo_level, n_groups, mean_purity,
                                mean_collision, mean_split, mean_coverage))

  list(per_level=per_level, summary=summary_tbl, h_corr=h_corr)
}


# ── PART B: PLOTTING ──────────────────────────────────────────────────────────

#' Lollipop / dot-range chart: purity, collision, split rate across pseudo levels
#'
#' @param pseudo_results  Output of run_pseudoclonotype_analysis().
#' @param subset_var      If not NULL, colour points by this column (e.g., "subset").
#' @return  ggplot2 object (patchwork of three panels).
plot_pseudoclonotype_summary <- function(pseudo_results, subset_var = NULL) {

  per <- pseudo_results$per_level
  smry <- pseudo_results$summary

  level_order <- levels(smry$pseudo_level)
  if (is.null(level_order))
    level_order <- unique(smry$pseudo_level)

  per <- per %>%
    mutate(pseudo_level = factor(pseudo_level, levels=level_order))
  smry <- smry %>%
    mutate(pseudo_level = factor(pseudo_level, levels=level_order))

  color_aes <- if (!is.null(subset_var) && subset_var %in% colnames(per))
    aes(color=.data[[subset_var]]) else aes()

  subset_scale <- if (!is.null(subset_var) && subset_var %in% colnames(per))
    scale_color_manual(values=c(CD8="#D6604D",CD4="#2166AC",
                                 Treg="#4DAF4A",DN="#FF7F00",Unknown="grey60"),
                        name="Subset")
  else NULL

  p1 <- ggplot(per, aes(x=pseudo_level, y=purity)) +
    geom_jitter(color_aes, width=.15, size=2.5, alpha=.7) +
    geom_point(data=smry, aes(y=mean_purity), color="black", size=4, shape=18) +
    geom_errorbar(data=smry, aes(y=mean_purity,
                                  ymin=mean_purity-sd_purity,
                                  ymax=mean_purity+sd_purity),
                  width=.2, linewidth=.8, color="black") +
    subset_scale +
    scale_y_continuous(limits=c(0,1), labels=percent_format()) +
    labs(title="Purity", subtitle="Higher = each pseudo maps to one real clone",
         x=NULL, y="Weighted purity") +
    theme_bw(base_size=11) +
    theme(axis.text.x=element_text(angle=35, hjust=1))

  p2 <- ggplot(per, aes(x=pseudo_level, y=collision_rate)) +
    geom_jitter(color_aes, width=.15, size=2.5, alpha=.7) +
    geom_point(data=smry, aes(y=mean_collision), color="black", size=4, shape=18) +
    subset_scale +
    scale_y_continuous(limits=c(0,1), labels=percent_format()) +
    labs(title="Collision Rate",
         subtitle="Fraction of pseudos containing >1 real clonotype",
         x=NULL, y="Collision rate") +
    theme_bw(base_size=11) +
    theme(axis.text.x=element_text(angle=35, hjust=1), legend.position="none")

  p3 <- ggplot(per, aes(x=pseudo_level, y=split_rate)) +
    geom_jitter(color_aes, width=.15, size=2.5, alpha=.7) +
    geom_point(data=smry, aes(y=mean_split), color="black", size=4, shape=18) +
    subset_scale +
    scale_y_continuous(limits=c(0,1), labels=percent_format()) +
    labs(title="Split Rate",
         subtitle="Fraction of real clonotypes spanning >1 pseudo",
         x=NULL, y="Split rate") +
    theme_bw(base_size=11) +
    theme(axis.text.x=element_text(angle=35, hjust=1), legend.position="none")

  (p1 | p2 | p3) +
    plot_annotation(
      title    = "Pseudoclonotype Recapitulation — Gene-Combo Level Comparison",
      subtitle = "Diamond = group mean ± SD. TRBV = coarsest; Full_paired = finest proxy.",
      theme    = theme(plot.title=element_text(face="bold",size=13))
    )
}

#' Scatter: H_pseudoclonotype vs H_real per group, faceted by pseudo level
#'
#' @param pseudo_results  Output of run_pseudoclonotype_analysis().
#' @param group_var       Column used to colour/label points.
#' @param subset_var      If not NULL, facet within each level by subset.
#' @return  ggplot2 object.
plot_pseudo_h_scatter <- function(pseudo_results,
                                   group_var  = "sample_id",
                                   subset_var = NULL) {

  per <- pseudo_results$per_level
  h_cor <- pseudo_results$h_corr

  level_order <- unique(per$pseudo_level)

  rho_labels <- h_cor %>%
    transmute(pseudo_level,
              label = paste0("ρ = ", sprintf("%.2f", spearman_rho)))

  per <- left_join(per, rho_labels, by="pseudo_level") %>%
    mutate(pseudo_level = factor(pseudo_level, levels=level_order))

  color_aes <- if (!is.null(subset_var) && subset_var %in% colnames(per))
    aes(color=.data[[subset_var]]) else aes(color=.data[[group_var]])

  all_h <- c(per$H_pseudo, per$H_real)
  h_rng <- range(all_h, na.rm=TRUE)

  ggplot(per, aes(x=H_pseudo, y=H_real)) +
    geom_abline(slope=1, intercept=0, color="red", linetype="dotted", linewidth=.7) +
    geom_point(color_aes, size=2.8, alpha=.8) +
    geom_smooth(method="lm", se=TRUE, color="grey40",
                linetype="dashed", linewidth=.7) +
    geom_text(aes(label=label), x=Inf, y=-Inf, hjust=1.05, vjust=-.3,
              size=3, color="grey30", inherit.aes=FALSE) +
    facet_wrap(vars(pseudo_level), scales="free", ncol=4) +
    labs(
      title    = "H(pseudoclonotype) vs H(real clonotype)",
      subtitle = "Each point = one sample × subset group. Red dotted = 1:1. ρ shown per facet.",
      x        = "Shannon diversity of pseudoclonotype labels",
      y        = "Shannon diversity of real CDR3 clonotypes"
    ) +
    theme_bw(base_size=11) +
    theme(strip.text=element_text(face="bold", size=9), legend.position="bottom")
}

#' Purity heatmap: pseudo level × subset (mean purity)
#'
#' @param pseudo_results  Output of run_pseudoclonotype_analysis().
#' @param subset_var      Column name of the subset variable.
#' @return  Invisible NULL (pheatmap renders directly).
plot_purity_heatmap <- function(pseudo_results, subset_var = "subset") {

  per <- pseudo_results$per_level

  if (!subset_var %in% colnames(per)) {
    warning("Column '", subset_var, "' not found; skipping purity heatmap.")
    return(invisible(NULL))
  }

  mat <- per %>%
    group_by(pseudo_level, .data[[subset_var]]) %>%
    summarise(mean_purity = mean(purity, na.rm=TRUE), .groups="drop") %>%
    pivot_wider(names_from=all_of(subset_var),
                values_from=mean_purity, values_fill=NA) %>%
    tibble::column_to_rownames("pseudo_level") %>%
    as.matrix()

  pheatmap::pheatmap(
    mat,
    color       = colorRampPalette(c("#fee2e2","#fef3c7","#d1fae5"))(100),
    main        = "Mean Pseudoclonotype Purity\n(by gene-combo level × T cell subset)",
    display_numbers = TRUE,
    number_format   = "%.2f",
    fontsize_row = 10, fontsize_col = 11,
    cluster_rows = FALSE, cluster_cols = FALSE,
    border_color = "white", na_col = "grey90"
  )
  invisible(NULL)
}


# =============================================================================
# SECTION 14: RNA-SEQ TCR GENE EXPRESSION ANALYSIS
# =============================================================================
#
# OVERVIEW:
#   Compares two sources of TCR gene usage information:
#
#   (A) TCR-seq calls  — V/J gene assignments from VDJ sequencing (precise,
#       sequence-level, one call per cell).
#   (B) RNA-seq expression — which TCR gene locus has the highest normalised
#       expression in the scRNA-seq matrix (noisier, but available in any
#       scRNA-seq dataset without matched VDJ library).
#
#   Concordance between A and B validates that RNA-seq expression can serve
#   as a fallback when TCR-seq is unavailable.
#
# SUBSET CLASSIFICATION (extended):
#   DP  — CD4 > 0  AND  CD8A/B > 0  (double positive)
#   CD8 — CD8A/B > 0  AND  NOT CD4  (single positive CD8)
#   Treg— CD4 > 0  AND  FOXP3 > 0  AND  NOT CD8
#   CD4 — CD4 > 0  AND  NOT CD8  AND  NOT FOXP3  (single positive CD4)
#   DN  — CD3D/E > 0  AND  NOT CD4  AND  NOT CD8  (double negative)
#   Non_T — no CD3 signal
#
# KEY FUNCTIONS:
#   load_h5ad_with_expression()
#   assign_rnaseq_dominant_gene()
#   compute_rnaseq_gene_usage_diversity()
#   compare_rnaseq_vs_tcrseq()
#   plot_tcr_expression_heatmap()
#   annotate_subsets_extended()        # adds DP and DN
# =============================================================================


# ── 14.1  Extended subset annotation (CD4 / CD8 / DP / DN / Treg) ────────────

# Extended annotation-string patterns
.SUBSET_MAP_EXT <- list(
  DP   = c("DP","double.positive","CD4.*CD8","CD8.*CD4","DP T","double pos"),
  CD8  = c("CD8","cytotoxic","CTL","cytotox","Tex","exhausted","effector CD8","CD8T"),
  Treg = c("Treg","T_reg","regulatory","FOXP3","T regulatory"),
  CD4  = c("CD4","T helper","Th1","Th2","Th17","Tfh","T_helper","CD4T","naive CD4"),
  DN   = c("double.negative","DN_T","DN T","gamma.delta","gdT","gd T","NKT")
)

#' Map raw annotation strings to extended canonical T cell subsets
#'
#' Returns one of: "DP", "CD8", "Treg", "CD4", "DN", "Non_T", "Unknown".
#' Order matters: DP is checked before CD8 and CD4.
#'
#' @param x  Character vector of raw annotations.
#' @return   Character vector of canonical subset labels.
standardise_subset_extended <- function(x) {
  x   <- as.character(x)
  out <- rep("Unknown", length(x))
  out[is.na(x) | x %in% c("","None","nan","NaN")] <- "Unknown"
  out[grepl("non.?t|nonT|non_T|B cell|NK[^T]|Mono|DC|macro|plasma|mast|neutro",
             x, ignore.case=TRUE)] <- "Non_T"
  for (label in c("DP","Treg","CD8","CD4","DN")) {
    pat <- paste(.SUBSET_MAP_EXT[[label]], collapse="|")
    out[grepl(pat, x, ignore.case=TRUE)] <- label
  }
  out
}

#' Classify cells using marker gene counts into DP/CD8/Treg/CD4/DN/Non_T
#'
#' Hierarchical rule applied to raw count columns (cd4_count, cd8a_count, etc.)
#' as produced by annotate_tcell_subsets_from_aggr().
#'
#' @param cell_dt   data.table with columns: cd4_count, cd8a_count, cd8b_count,
#'                  cd3d_count, cd3e_count, foxp3_count.
#' @param thr       Minimum raw UMI count to call a marker positive (default 1).
#' @return  cell_dt with a 'subset' column added/replaced.
classify_subsets_extended <- function(cell_dt, thr = 1) {
  cell_dt[, subset := fcase(
    # DP: CD4 AND (CD8A OR CD8B)
    (cd4_count >= thr) & (cd8a_count >= thr | cd8b_count >= thr),      "DP",
    # CD8: CD8 only
    (cd8a_count >= thr | cd8b_count >= thr),                            "CD8",
    # Treg: CD4 + FOXP3, no CD8
    (cd4_count >= thr) & (foxp3_count >= thr),                          "Treg",
    # CD4: CD4 only
    (cd4_count >= thr),                                                  "CD4",
    # DN: CD3 signal but no CD4 or CD8
    (cd3d_count >= thr | cd3e_count >= thr),                            "DN",
    default = "Non_T"
  )]
  cell_dt[]
}

#' Add extended subset annotation (DP / CD8 / Treg / CD4 / DN) to TCR metadata
#'
#' Replaces add_subset_from_annotation() with extended classification.
#' Detects annotation column and maps to the 6-way scheme.
#'
#' @param metadata     Data frame with cell-level metadata.
#' @param subset_col   Override: exact column to use (NULL = auto-detect).
#' @param keep_nont    If FALSE (default), drop Non_T cells.
#' @param keep_unknown If FALSE (default), drop Unknown cells.
#' @return  metadata with 'subset' column (one of DP/CD8/Treg/CD4/DN/Non_T/Unknown).
add_subset_extended <- function(metadata,
                                 subset_col   = NULL,
                                 keep_nont    = FALSE,
                                 keep_unknown = FALSE) {

  if (is.null(subset_col)) subset_col <- detect_subset_column(metadata)

  if (is.null(subset_col)) {
    warning("No subset annotation column detected. All cells labelled 'Unknown'.")
    metadata$subset <- "Unknown"
    return(metadata)
  }

  message("  Using '", subset_col, "' for extended T cell subset annotation.")
  metadata$subset <- standardise_subset_extended(metadata[[subset_col]])

  tab <- table(metadata$subset)
  message("  Subset mapping (DP/CD8/Treg/CD4/DN):")
  for (nm in names(tab))
    message(sprintf("    %-10s: %d (%.1f%%)", nm, tab[nm], 100*tab[nm]/nrow(metadata)))

  if (!keep_nont)    metadata <- metadata[metadata$subset != "Non_T",    ]
  if (!keep_unknown) metadata <- metadata[metadata$subset != "Unknown",  ]
  metadata
}


# ── 14.2  Load h5ad with expression matrix ────────────────────────────────────

#' Load an .h5ad file returning both cell metadata and the expression matrix
#'
#' The expression matrix is filtered to only the genes in `genes_keep` before
#' returning, to avoid loading the full matrix into memory unnecessarily.
#'
#' @param path        Path to the .h5ad file.
#' @param genes_keep  Character vector of gene symbols to retain (e.g.,
#'                    c(TCR_GENES$TRBV, TCR_GENES$TRAV, ...)).
#'                    If NULL, the full matrix is returned.
#' @param layer       Assay/layer name to use. NULL tries logcounts → counts →
#'                    first available assay (zellkonverter) or X (anndata).
#' @return  List with:
#'   $meta  — data frame (one row per cell, barcode_aggr as first column)
#'   $expr  — sparse or dense matrix (genes × cells), rownames = gene symbols
#'   $genes — character vector of genes actually found in the matrix
load_h5ad_with_expression <- function(path,
                                       genes_keep = NULL,
                                       layer      = NULL) {
  stopifnot(file.exists(path))
  message("  Loading expression matrix from: ", basename(path))

  if (requireNamespace("zellkonverter", quietly=TRUE)) {
    sce  <- zellkonverter::readH5AD(path, X=TRUE, use_hdf5=FALSE)
    meta <- as.data.frame(SingleCellExperiment::colData(sce))
    meta <- tibble::rownames_to_column(meta, "barcode_aggr")

    # Pick assay
    avail <- SummarizedExperiment::assayNames(sce)
    chosen <- if (!is.null(layer) && layer %in% avail) layer
              else if ("logcounts" %in% avail)          "logcounts"
              else if ("lognorm"   %in% avail)          "lognorm"
              else if ("counts"    %in% avail)          "counts"
              else                                       avail[1]
    message("  Using assay: '", chosen, "'")
    expr <- SummarizedExperiment::assay(sce, chosen)

    found <- character(0)
    if (!is.null(genes_keep)) {
      keep  <- intersect(genes_keep, rownames(expr))
      expr  <- expr[keep, , drop=FALSE]
      found <- keep
      message("  TCR genes found in matrix: ", length(keep),
              " / ", length(genes_keep))
    } else {
      found <- rownames(expr)
    }
    return(list(meta=meta, expr=expr, genes=found))
  }

  if (requireNamespace("anndata", quietly=TRUE)) {
    ad   <- anndata::read_h5ad(path)
    meta <- as.data.frame(ad$obs)
    meta <- tibble::rownames_to_column(meta, "barcode_aggr")
    expr <- t(ad$X)           # genes × cells
    rownames(expr) <- ad$var_names

    found <- character(0)
    if (!is.null(genes_keep)) {
      keep  <- intersect(genes_keep, rownames(expr))
      expr  <- expr[keep, , drop=FALSE]
      found <- keep
    } else {
      found <- rownames(expr)
    }
    return(list(meta=meta, expr=expr, genes=found))
  }

  stop("Install zellkonverter (BiocManager) or anndata (CRAN).")
}


# ── 14.3  RNA-seq gene call assignment ────────────────────────────────────────

#' Assign dominant TCR gene per cell from scRNA-seq expression
#'
#' For each cell (column), finds the canonical TCR gene (from gene_list) with
#' the highest expression value. Cells with no detected TCR gene expression
#' (max value <= min_expr) are returned as NA.
#'
#' @param expr       Matrix (genes × cells) of normalised expression values.
#'                   Rownames must be gene symbols.
#' @param gene_list  Character vector of canonical genes to consider
#'                   (e.g., TCR_GENES$TRBV). Genes absent from expr are skipped.
#' @param min_expr   Minimum expression to call a gene as detected (default 0,
#'                   i.e., any positive value).
#' @return  Named character vector (length = ncol(expr)); names = cell barcodes.
assign_rnaseq_dominant_gene <- function(expr, gene_list, min_expr = 0) {
  keep <- intersect(gene_list, rownames(expr))
  if (length(keep) == 0) {
    warning("No genes from gene_list found in expression matrix rownames.")
    return(setNames(rep(NA_character_, ncol(expr)), colnames(expr)))
  }

  sub <- expr[keep, , drop=FALSE]

  # Handle sparse matrices (Matrix package)
  if (inherits(sub, "sparseMatrix"))
    sub <- as.matrix(sub)

  max_vals <- apply(sub, 2, max)
  max_idx  <- apply(sub, 2, which.max)

  calls <- keep[max_idx]
  calls[max_vals <= min_expr] <- NA_character_
  setNames(calls, colnames(expr))
}

#' Add RNA-seq-derived TCR gene calls to cell metadata
#'
#' Assigns a dominant TRBV gene (and optionally TRAV, TRBJ, TRAJ) per cell
#' from the expression matrix and appends as new columns:
#'   trbv_rna, trav_rna, trbj_rna, traj_rna
#'
#' @param meta       Data frame with 'barcode_aggr' column.
#' @param expr_data  Output of load_h5ad_with_expression() ($expr and $genes).
#' @param min_expr   Minimum expression value to call a gene detected.
#' @return  meta with rna_* columns added.
add_rnaseq_gene_calls <- function(meta, expr_data, min_expr = 0) {
  expr  <- expr_data$expr

  # Make sure barcodes align
  bc    <- meta$barcode_aggr
  bc_in <- intersect(bc, colnames(expr))
  if (length(bc_in) == 0) {
    warning("No barcode overlap between metadata and expression matrix.")
    for (col in c("trbv_rna","trav_rna","trbj_rna","traj_rna"))
      meta[[col]] <- NA_character_
    return(meta)
  }
  message("  Barcode overlap: ", length(bc_in), " / ", length(bc), " cells")

  chain_defs <- list(
    trbv_rna = TCR_GENES$TRBV,
    trav_rna = TCR_GENES$TRAV,
    trbj_rna = TCR_GENES$TRBJ,
    traj_rna = TCR_GENES$TRAJ
  )

  for (col_name in names(chain_defs)) {
    gene_list <- chain_defs[[col_name]]
    calls     <- assign_rnaseq_dominant_gene(
                   expr[, bc_in, drop=FALSE], gene_list, min_expr)
    call_vec  <- setNames(rep(NA_character_, nrow(meta)), meta$barcode_aggr)
    call_vec[bc_in] <- calls[bc_in]
    meta[[col_name]] <- unname(call_vec[meta$barcode_aggr])
    n_called <- sum(!is.na(meta[[col_name]]))
    message(sprintf("  %-12s: %d / %d cells assigned (%.1f%%)",
                    col_name, n_called, nrow(meta),
                    100*n_called/nrow(meta)))
  }
  meta
}


# ── 14.4  RNA-seq vs TCR-seq comparison ──────────────────────────────────────

#' Compare RNA-seq and TCR-seq gene usage per sample × subset
#'
#' Computes:
#'   concordance      — fraction of cells where RNA-seq call matches TCR-seq call
#'   coverage_rna     — fraction of cells with RNA-seq call
#'   coverage_tcrseq  — fraction of cells with TCR-seq call
#'   H_rna / H_tcrseq — Shannon diversity of usage (per group)
#'   spearman_rho     — correlation between H_rna and H_tcrseq across groups
#'
#' @param metadata   Data frame with columns trbv_gene (TCR-seq) and trbv_rna
#'                   (RNA-seq), plus group_vars.
#' @param rna_col    RNA-seq gene call column name (default "trbv_rna").
#' @param tcr_col    TCR-seq gene call column name (default "trbv_gene").
#' @param group_vars Grouping columns (e.g., c("sample_id","subset")).
#' @param min_cells  Minimum cells per group.
#' @return  List with $per_group (data frame) and $correlations (data frame).
compare_rnaseq_vs_tcrseq <- function(metadata,
                                      rna_col    = "trbv_rna",
                                      tcr_col    = "trbv_gene",
                                      group_vars = "sample_id",
                                      min_cells  = 30) {
  for (col in c(rna_col, tcr_col))
    if (!col %in% colnames(metadata))
      stop("Column '", col, "' not found in metadata.")

  per_group <- metadata %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      n_cells         = n(),
      n_with_rna      = sum(!is.na(.data[[rna_col]]) & .data[[rna_col]] != ""),
      n_with_tcr      = sum(!is.na(.data[[tcr_col]]) & .data[[tcr_col]] != ""),
      n_both          = sum(!is.na(.data[[rna_col]]) & .data[[rna_col]] != "" &
                             !is.na(.data[[tcr_col]]) & .data[[tcr_col]] != ""),
      coverage_rna    = round(n_with_rna / n_cells, 3),
      coverage_tcrseq = round(n_with_tcr / n_cells, 3),
      # Concordance: only among cells with both calls
      concordance     = {
        both <- !is.na(.data[[rna_col]]) & .data[[rna_col]] != "" &
                !is.na(.data[[tcr_col]]) & .data[[tcr_col]] != ""
        if (sum(both) == 0) NA_real_
        else round(mean(.data[[rna_col]][both] == .data[[tcr_col]][both]), 3)
      },
      H_rna    = {
        vals <- .data[[rna_col]][!is.na(.data[[rna_col]]) & .data[[rna_col]] != ""]
        if (length(vals) < 2) NA_real_
        else shannon_diversity(as.numeric(table(vals)))
      },
      H_tcrseq = {
        vals <- .data[[tcr_col]][!is.na(.data[[tcr_col]]) & .data[[tcr_col]] != ""]
        if (length(vals) < 2) NA_real_
        else shannon_diversity(as.numeric(table(vals)))
      },
      .groups = "drop"
    ) %>%
    filter(n_cells >= min_cells) %>%
    mutate(H_diff = H_rna - H_tcrseq)

  # Correlations
  correlations <- list()
  ok <- !is.na(per_group$H_rna) & !is.na(per_group$H_tcrseq)
  if (sum(ok) >= 3) {
    correlations$spearman <- cor.test(per_group$H_rna[ok],
                                       per_group$H_tcrseq[ok],
                                       method="spearman", exact=FALSE)
    correlations$pearson  <- cor.test(per_group$H_rna[ok],
                                       per_group$H_tcrseq[ok],
                                       method="pearson")
    message(sprintf(
      "  H_RNA vs H_TCRseq: Spearman rho=%.3f (p=%.4f), n=%d",
      correlations$spearman$estimate, correlations$spearman$p.value, sum(ok)
    ))
    message(sprintf("  Mean concordance: %.1f%%",
                    100*mean(per_group$concordance, na.rm=TRUE)))
  }

  # Per-gene concordance
  shared <- metadata %>%
    filter(!is.na(.data[[rna_col]]), .data[[rna_col]] != "",
           !is.na(.data[[tcr_col]]), .data[[tcr_col]] != "")

  if (nrow(shared) > 0) {
    gene_concordance <- shared %>%
      group_by(tcr_gene = .data[[tcr_col]]) %>%
      summarise(
        n_cells       = n(),
        n_concordant  = sum(.data[[rna_col]] == .data[[tcr_col]]),
        concordance   = round(n_concordant / n_cells, 3),
        top_rna_call  = names(sort(table(.data[[rna_col]]), decreasing=TRUE))[1],
        .groups = "drop"
      ) %>%
      arrange(desc(n_cells))
  } else {
    gene_concordance <- tibble()
  }

  list(per_group=per_group, correlations=correlations,
       gene_concordance=gene_concordance)
}


# ── 14.5  TCR gene expression heatmap ─────────────────────────────────────────

#' Heatmap of average TCR gene expression per sample × subset group
#'
#' Computes mean normalised expression of each detected TCR gene across cells in
#' each group, then plots a clustered heatmap (genes × groups).
#'
#' @param expr        Matrix (genes × cells) from load_h5ad_with_expression().
#' @param meta        Data frame with barcode_aggr and group_vars columns.
#' @param group_vars  Grouping columns (e.g., c("sample_id","subset")).
#' @param chains      Character vector of chain prefixes to include.
#'                    Default c("TRBV","TRAV","TRBJ","TRAJ").
#' @param min_pct     Minimum fraction of cells with detectable expression
#'                    for a gene to be shown (default 0.01 = 1%).
#' @param title       Plot title.
#' @return  Invisible NULL (pheatmap renders directly); also returns the
#'          numeric matrix invisibly for downstream use.
plot_tcr_expression_heatmap <- function(expr,
                                         meta,
                                         group_vars = c("sample_id","subset"),
                                         chains     = c("TRBV","TRAV","TRBJ","TRAJ"),
                                         min_pct    = 0.01,
                                         title      = "TCR Gene Expression\n(mean log-normalised per group)") {

  # Keep TCR genes present in expression matrix
  tcr_all  <- unlist(TCR_GENES[intersect(chains, names(TCR_GENES))], use.names=FALSE)
  tcr_keep <- intersect(tcr_all, rownames(expr))

  if (length(tcr_keep) == 0) {
    warning("No canonical TCR genes found in expression matrix rownames.")
    return(invisible(NULL))
  }
  message("  TCR genes in matrix: ", length(tcr_keep))

  sub_expr <- expr[tcr_keep, , drop=FALSE]

  # Filter genes expressed in at least min_pct of cells
  pct_expr <- Matrix::rowMeans(sub_expr > 0)
  sub_expr <- sub_expr[pct_expr >= min_pct, , drop=FALSE]
  if (nrow(sub_expr) == 0) {
    warning("No TCR genes pass min_pct filter (", min_pct, ").")
    return(invisible(NULL))
  }
  message("  Genes passing min_pct (", min_pct, "): ", nrow(sub_expr))

  # Align barcodes
  bc_common <- intersect(meta$barcode_aggr, colnames(sub_expr))
  if (length(bc_common) == 0) {
    warning("No barcode overlap between metadata and expression matrix.")
    return(invisible(NULL))
  }
  meta_sub  <- meta[meta$barcode_aggr %in% bc_common, , drop=FALSE]
  expr_sub  <- sub_expr[, bc_common, drop=FALSE]

  # Build group labels
  meta_sub$group_label <- apply(
    meta_sub[, group_vars, drop=FALSE], 1,
    function(x) paste(x, collapse=" | ")
  )

  # Compute mean expression per group
  groups <- unique(meta_sub$group_label)
  mean_mat <- sapply(groups, function(g) {
    bc_g <- meta_sub$barcode_aggr[meta_sub$group_label == g]
    if (length(bc_g) == 0) return(rep(NA_real_, nrow(expr_sub)))
    rowMeans(as.matrix(expr_sub[, bc_g, drop=FALSE]), na.rm=TRUE)
  })
  rownames(mean_mat) <- rownames(expr_sub)
  colnames(mean_mat) <- groups

  # Drop all-NA columns
  mean_mat <- mean_mat[, colSums(!is.na(mean_mat)) > 0, drop=FALSE]

  # Annotation: colour bars by chain
  chain_map <- sapply(rownames(mean_mat), function(g) {
    ch <- sub("^(TRBV|TRAV|TRBJ|TRAJ|TRBD).*","\\1", g)
    if (ch %in% c("TRBV","TRAV","TRBJ","TRAJ","TRBD")) ch else "Other"
  })
  row_ann <- data.frame(Chain=chain_map, row.names=rownames(mean_mat))
  ann_col <- list(Chain=c(TRBV="#2166AC",TRAV="#EC4899",
                            TRBJ="#22C55E",TRAJ="#A855F7",TRBD="#F59E0B"))

  pheatmap::pheatmap(
    mean_mat,
    color            = colorRampPalette(c("white","#FEF3C7","#F59E0B","#DC2626"))(100),
    main             = title,
    annotation_row   = row_ann,
    annotation_colors= ann_col,
    cluster_rows     = TRUE,
    cluster_cols     = TRUE,
    fontsize_row     = 6,
    fontsize_col     = 8,
    border_color     = NA,
    na_col           = "grey90",
    show_rownames    = TRUE
  )
  invisible(mean_mat)
}


# ── 14.6  Scatter plots for RNA-seq vs TCR-seq comparison ────────────────────

#' Scatter: H_RNA vs H_TCRseq per group, coloured by subset or sample
#'
#' @param comparison_result  Output of compare_rnaseq_vs_tcrseq().
#' @param group_var          Column used to colour/label points.
#' @param color_var          Column used for colour (default = group_var).
#' @param chain_label        Label for the chain analysed (e.g., "TRBV").
#' @return  ggplot2 object.
plot_rnaseq_vs_tcrseq_scatter <- function(comparison_result,
                                           group_var   = "sample_id",
                                           color_var   = NULL,
                                           chain_label = "TRBV") {
  pg <- comparison_result$per_group
  if (is.null(color_var)) color_var <- group_var

  sp_txt <- ""
  if (!is.null(comparison_result$correlations$spearman)) {
    sp  <- comparison_result$correlations$spearman
    sp_txt <- paste0("Spearman ρ = ", round(sp$estimate,3),
                     "  (p = ", signif(sp$p.value,3), ")")
  }

  ggplot(pg, aes(H_tcrseq, H_rna, color=.data[[color_var]], label=.data[[group_var]])) +
    geom_abline(slope=1, intercept=0, color="red", linetype="dotted", linewidth=.7) +
    geom_smooth(method="lm", se=TRUE, color="grey40",
                linetype="dashed", linewidth=.7) +
    geom_point(size=3, alpha=.85) +
    geom_text_repel(size=2.8, max.overlaps=15) +
    scale_color_manual(
      values=c(CD8="#D6604D",CD4="#2166AC",Treg="#4DAF4A",DP="#FF7F00",
                DN="#9B59B6",Unknown="grey60"),
      na.value="grey60"
    ) +
    labs(
      title    = paste(chain_label, "Usage Diversity: RNA-seq vs TCR-seq"),
      subtitle = sp_txt,
      x        = paste0("H_", chain_label, " from TCR-seq V-gene calls"),
      y        = paste0("H_", chain_label, " from RNA-seq expression"),
      caption  = "Red dotted = 1:1. Points above = RNA-seq overestimates diversity."
    ) +
    theme_bw(base_size=12) + theme(legend.position="bottom")
}

#' Bar chart: per-gene concordance between RNA-seq and TCR-seq calls
#'
#' Shows what fraction of cells assigned to each TCR-seq gene have the same
#' gene as their top RNA-seq expressed gene.
#'
#' @param comparison_result  Output of compare_rnaseq_vs_tcrseq().
#' @param top_n              Number of top genes (by cell count) to show.
#' @param chain_label        Label for the chain (e.g., "TRBV").
#' @return  ggplot2 object.
plot_gene_concordance_bar <- function(comparison_result,
                                       top_n       = 20,
                                       chain_label = "TRBV") {
  gc <- comparison_result$gene_concordance

  if (nrow(gc) == 0) {
    warning("No gene concordance data available.")
    return(ggplot() + annotate("text",x=.5,y=.5,label="No data") + theme_void())
  }

  gc_plot <- gc %>%
    slice_max(n_cells, n=top_n) %>%
    mutate(tcr_gene = factor(tcr_gene, levels=rev(tcr_gene)))

  ggplot(gc_plot, aes(x=concordance, y=tcr_gene, fill=concordance)) +
    geom_col(width=.7) +
    geom_text(aes(label=paste0(round(concordance*100,0),"%")),
              hjust=-.05, size=3) +
    scale_fill_gradient2(low="#fee2e2", mid="#fef3c7", high="#d1fae5",
                          midpoint=.5, limits=c(0,1), guide="none") +
    scale_x_continuous(limits=c(0,1.1), labels=percent_format()) +
    labs(
      title    = paste(chain_label, "Per-Gene Concordance: RNA-seq vs TCR-seq"),
      subtitle = "Fraction of cells where RNA-seq dominant gene matches TCR-seq call",
      x        = "Concordance rate", y = chain_label
    ) +
    theme_bw(base_size=11)
}
