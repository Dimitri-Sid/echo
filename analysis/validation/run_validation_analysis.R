# =============================================================================
# Validation Analysis — TRBV Proxy vs True Clonotype Diversity
# Published Datasets: Blum, Lechner, Luoma
# =============================================================================
#
# DATASETS (all stored as .h5ad AnnData objects):
#   blum_combined_tissue_data.h5ad
#     Blum et al. — combined tissue T cells (multiple tissues/donors)
#   lechner_thyroiditis_tissue_cd8_seurat_object.h5ad
#     Lechner et al. — thyroiditis tissue CD8 T cells
#   luoma_colitis_cd8_seurat_object.h5ad
#     Luoma et al. — colitis tissue CD8 T cells
#
# WHAT THIS SCRIPT DOES:
#   1. Loads each .h5ad file via zellkonverter (or anndata as fallback).
#   2. Auto-detects and normalises TCR column names to the schema expected
#      by tcr_diversity_functions.R (trbv_gene, clonotype_id, sample_id).
#   3. Applies the full NanoString-style TCR diversity pipeline to each dataset.
#   4. Generates all comparison plots and summary tables.
#   5. Produces a cross-dataset summary comparing how well the TRBV proxy
#      tracks true clonotype diversity across contexts.
#
# OUTPUTS (written to output/validation/):
#   Per dataset:
#     <dataset>_tcr_diversity_summary.csv
#     <dataset>_trbv_stacked_bar.pdf
#     <dataset>_clonotype_rank_abundance.pdf
#     <dataset>_h_trbv_vs_h_clonotype.pdf
#     <dataset>_clonality_comparison.pdf
#     <dataset>_trbv_heatmap.pdf
#     <dataset>_within_trbv_diversity.pdf
#     <dataset>_bland_altman.pdf
#     <dataset>_failure_mode_flags.csv
#   Cross-dataset:
#     cross_dataset_correlation_summary.csv
#     cross_dataset_h_trbv_vs_h_clonotype.pdf
# =============================================================================


# ── 0. Setup ------------------------------------------------------------------

# Install missing packages (run once):
# BiocManager::install("zellkonverter")
# install.packages(c("dplyr","tidyr","purrr","stringr","ggplot2","ggrepel",
#                    "scales","data.table","RColorBrewer","patchwork","pheatmap"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(data.table)
  library(RColorBrewer)
  library(patchwork)
  library(pheatmap)
})

REPO_ROOT   <- here::here()
source(file.path(REPO_ROOT, "R", "tcr_diversity_functions.R"))

VAL_DIR    <- file.path(REPO_ROOT, "validation")
OUTPUT_DIR <- file.path(REPO_ROOT, "output", "validation")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

MIN_CELLS    <- 30
DOWNSAMPLE_N <- 200   # set to 0 to skip sensitivity analysis


# ── 1. h5ad loading helpers ---------------------------------------------------

#' Load an .h5ad file and return the cell metadata as a data frame.
#'
#' Tries zellkonverter first (Bioconductor), then falls back to the anndata
#' R package (CRAN). Returns a data frame with one row per cell; rownames
#' are cell barcodes (moved to column "barcode_aggr").
load_h5ad_metadata <- function(path) {
  stopifnot(file.exists(path))
  message("  Loading: ", basename(path))

  # --- zellkonverter path (preferred) ---
  if (requireNamespace("zellkonverter", quietly = TRUE)) {
    sce  <- zellkonverter::readH5AD(path, X = FALSE, use_hdf5 = FALSE)
    meta <- as.data.frame(SingleCellExperiment::colData(sce))
    meta <- tibble::rownames_to_column(meta, "barcode_aggr")
    message("  -> ", nrow(meta), " cells loaded via zellkonverter")
    return(meta)
  }

  # --- anndata path (fallback) ---
  if (requireNamespace("anndata", quietly = TRUE)) {
    ad   <- anndata::read_h5ad(path)
    meta <- as.data.frame(ad$obs)
    meta <- tibble::rownames_to_column(meta, "barcode_aggr")
    message("  -> ", nrow(meta), " cells loaded via anndata")
    return(meta)
  }

  stop(
    "Neither 'zellkonverter' (Bioconductor) nor 'anndata' (CRAN) is installed.\n",
    "Install one with:\n",
    "  BiocManager::install('zellkonverter')   # preferred\n",
    "  install.packages('anndata')             # fallback"
  )
}


# ── 2. TCR column auto-detection and normalisation ----------------------------

# Known column-name patterns (regex, case-insensitive) for each required field.
# First match wins.
.TRBV_PATTERNS <- c(
  "^v_gene$", "^trbv$", "^trbv_gene$", "^v_gene_beta$",
  "IR_VDJ_1_v_call", "^TCR_Beta_V_gene$", "^beta_v_gene$",
  "^vb_gene$", "^TRBV_gene$", "v_call"
)

.CLONO_PATTERNS <- c(
  "^clonotype_id$", "^clone_id$", "^cloneid$",
  "^cdr3_beta$", "^cdr3_b_aa$", "^IR_obs_vj_1_junction_aa$",
  "^TCR_clonotype$", "^clonotype$", "^clone$",
  "^raw_clonotype_id$", "^ct_id$"
)

.SAMPLE_PATTERNS <- c(
  "^sample_id$", "^sample$", "^donor$", "^patient$",
  "^patient_id$", "^donor_id$", "^orig.ident$",
  "^batch$", "^library_id$"
)

.find_col <- function(colnames_vec, patterns) {
  for (pat in patterns) {
    hit <- grep(pat, colnames_vec, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  NULL
}

#' Auto-detect TCR columns and return a normalised metadata data frame.
#'
#' Adds/renames columns to: trbv_gene, clonotype_id, sample_id.
#' Prints a column-mapping report so you can verify the detection.
#'
#' @param meta       Data frame from load_h5ad_metadata().
#' @param dataset_id Short name used as fallback sample_id if no sample column exists.
#' @param trbv_col   Override: exact column name for TRBV gene (skips auto-detect).
#' @param clono_col  Override: exact column name for clonotype ID.
#' @param sample_col Override: exact column name for sample ID.
#' @return  Normalised data frame with trbv_gene, clonotype_id, sample_id.
normalise_tcr_metadata <- function(meta,
                                   dataset_id,
                                   trbv_col   = NULL,
                                   clono_col  = NULL,
                                   sample_col = NULL) {

  cols <- colnames(meta)

  # Auto-detect
  trbv_col   <- trbv_col   %||% .find_col(cols, .TRBV_PATTERNS)
  clono_col  <- clono_col  %||% .find_col(cols, .CLONO_PATTERNS)
  sample_col <- sample_col %||% .find_col(cols, .SAMPLE_PATTERNS)

  cat(sprintf(
    "\n  [%s] Column mapping:\n    TRBV gene    -> %s\n    Clonotype ID -> %s\n    Sample ID    -> %s\n",
    dataset_id,
    trbv_col  %||% "NOT FOUND",
    clono_col %||% "NOT FOUND",
    sample_col %||% paste0("NOT FOUND (using '", dataset_id, "')")
  ))

  if (is.null(trbv_col)) {
    warning("[", dataset_id, "] No TRBV gene column detected. ",
            "Set trbv_col= manually.\nAvailable columns: ",
            paste(cols, collapse = ", "))
  }
  if (is.null(clono_col)) {
    warning("[", dataset_id, "] No clonotype_id column detected. ",
            "Set clono_col= manually.\nAvailable columns: ",
            paste(cols, collapse = ", "))
  }

  # Rename to standard schema
  if (!is.null(trbv_col)  && trbv_col  != "trbv_gene")
    meta <- meta %>% rename(trbv_gene   = all_of(trbv_col))
  if (!is.null(clono_col) && clono_col != "clonotype_id")
    meta <- meta %>% rename(clonotype_id = all_of(clono_col))

  # sample_id: rename existing column or create constant
  if (!is.null(sample_col) && sample_col != "sample_id") {
    meta <- meta %>% rename(sample_id = all_of(sample_col))
  } else if (is.null(sample_col)) {
    meta$sample_id <- dataset_id
    message("  -> No sample column found; all cells labelled as '", dataset_id, "'")
  }

  # Ensure required columns exist (with NA if truly absent)
  if (!"trbv_gene"    %in% colnames(meta)) meta$trbv_gene    <- NA_character_
  if (!"clonotype_id" %in% colnames(meta)) meta$clonotype_id <- NA_character_

  # Tag with dataset name for cross-dataset plots
  meta$dataset <- dataset_id

  meta
}

# Null-coalescing operator
`%||%` <- function(a, b) if (!is.null(a)) a else b


# ── 3. Per-dataset analysis wrapper ------------------------------------------

#' Load, normalise, and run the full TCR diversity pipeline on one dataset.
#'
#' @param h5ad_path    Path to the .h5ad file.
#' @param dataset_id   Short name for outputs (e.g., "blum").
#' @param group_vars   Grouping variables for the analysis (default "sample_id").
#' @param trbv_col     Optional override for TRBV column name.
#' @param clono_col    Optional override for clonotype column name.
#' @param sample_col   Optional override for sample ID column name.
#' @param min_cells    Minimum cells per group.
#' @param downsample_n Cells to downsample to per group (0 = skip).
#' @return  List with $meta (normalised metadata), $results (diversity results),
#'          $summary_table, $benchmark.
run_dataset_analysis <- function(h5ad_path,
                                  dataset_id,
                                  group_vars   = "sample_id",
                                  trbv_col     = NULL,
                                  clono_col    = NULL,
                                  sample_col   = NULL,
                                  min_cells    = MIN_CELLS,
                                  downsample_n = DOWNSAMPLE_N) {

  message("\n", strrep("=", 70))
  message("DATASET: ", dataset_id)
  message(strrep("=", 70))

  # Load and inspect
  meta_raw <- load_h5ad_metadata(h5ad_path)

  cat("\nAll available columns:\n")
  cat(paste(sprintf("  [%03d] %s", seq_along(colnames(meta_raw)), colnames(meta_raw)),
            collapse = "\n"), "\n")

  # Normalise column names
  meta <- normalise_tcr_metadata(meta_raw, dataset_id,
                                  trbv_col, clono_col, sample_col)

  # Report coverage
  n_total  <- nrow(meta)
  n_trbv   <- sum(!is.na(meta$trbv_gene) & meta$trbv_gene != "" &
                    meta$trbv_gene != "None")
  n_clono  <- sum(!is.na(meta$clonotype_id) & meta$clonotype_id != "" &
                    meta$clonotype_id != "None")

  cat(sprintf("\n  Cells total:              %d\n", n_total))
  cat(sprintf("  Cells with TRBV gene:     %d (%.1f%%)\n",
              n_trbv, 100 * n_trbv / n_total))
  cat(sprintf("  Cells with clonotype ID:  %d (%.1f%%)\n",
              n_clono, 100 * n_clono / n_total))
  cat(sprintf("  Unique sample_ids:        %d\n",
              n_distinct(meta$sample_id)))
  cat(sprintf("  Unique TRBV genes:        %d\n",
              n_distinct(meta$trbv_gene[!is.na(meta$trbv_gene)], na.rm = TRUE)))
  cat(sprintf("  Unique clonotypes:        %d\n\n",
              n_distinct(meta$clonotype_id[!is.na(meta$clonotype_id)], na.rm = TRUE)))

  # Run diversity analysis
  results <- run_tcr_diversity_analysis(
    metadata     = meta,
    group_vars   = group_vars,
    min_cells    = min_cells,
    downsample_n = downsample_n
  )

  print_correlation_summary(results$comparison)

  summary_tbl <- build_summary_table(results$comparison)
  cat("\n--- Summary table ---\n")
  print(summary_tbl)

  # Save outputs
  out_prefix <- file.path(OUTPUT_DIR, dataset_id)

  write.csv(summary_tbl,
            paste0(out_prefix, "_tcr_diversity_summary.csv"),
            row.names = FALSE)

  benchmark <- results$benchmark
  if ("failure_mode" %in% colnames(benchmark$failures)) {
    write.csv(benchmark$failures,
              paste0(out_prefix, "_failure_mode_flags.csv"),
              row.names = FALSE)
  }
  write.csv(benchmark$within_trbv,
            paste0(out_prefix, "_within_trbv_clonotype_details.csv"),
            row.names = FALSE)

  # Plots
  message("  Generating plots for ", dataset_id, "...")

  group_var1 <- group_vars[1]   # primary axis for plots

  # a. TRBV stacked bar
  tryCatch({
    p <- plot_trbv_stacked_bar(meta, group_var = group_var1, top_n_genes = 20)
    ggsave(paste0(out_prefix, "_trbv_stacked_bar.pdf"), p, width = 10, height = 6)
  }, error = function(e) warning("trbv_stacked_bar failed: ", e$message))

  # b. Clonotype rank-abundance
  tryCatch({
    p <- plot_clonotype_rank_abundance(meta, group_var = group_var1)
    ggsave(paste0(out_prefix, "_clonotype_rank_abundance.pdf"), p, width = 9, height = 6)
  }, error = function(e) warning("rank_abundance failed: ", e$message))

  # c. H_TRBV vs H_clonotype scatter
  tryCatch({
    p <- plot_h_trbv_vs_h_clonotype(results$comparison$table,
                                     group_var = group_var1)
    ggsave(paste0(out_prefix, "_h_trbv_vs_h_clonotype.pdf"), p, width = 7, height = 6)
  }, error = function(e) warning("h_trbv_vs_h_clonotype failed: ", e$message))

  # d. Clonality comparison
  tryCatch({
    p <- plot_clonality_comparison(results$comparison$table,
                                    group_var = group_var1)
    ggsave(paste0(out_prefix, "_clonality_comparison.pdf"), p, width = 7, height = 6)
  }, error = function(e) warning("clonality_comparison failed: ", e$message))

  # e. TRBV heatmap
  tryCatch({
    pdf(paste0(out_prefix, "_trbv_heatmap.pdf"), width = 9, height = 8)
    plot_trbv_heatmap(meta, group_var = group_var1, top_n_genes = 25)
    dev.off()
  }, error = function(e) { dev.off(); warning("trbv_heatmap failed: ", e$message) })

  # f. Within-TRBV diversity (failure mode)
  tryCatch({
    p <- plot_within_trbv_diversity(benchmark, group_var = group_var1)
    ggsave(paste0(out_prefix, "_within_trbv_diversity.pdf"), p, width = 14, height = 8)
  }, error = function(e) warning("within_trbv_diversity failed: ", e$message))

  # g. Bland-Altman
  tryCatch({
    p <- plot_bland_altman(results$comparison$table, group_var = group_var1)
    ggsave(paste0(out_prefix, "_bland_altman.pdf"), p, width = 7, height = 6)
  }, error = function(e) warning("bland_altman failed: ", e$message))

  message("  Plots saved with prefix: ", out_prefix)

  list(
    meta          = meta,
    results       = results,
    summary_table = summary_tbl,
    benchmark     = benchmark
  )
}


# ── 4. Run all three validation datasets -------------------------------------

# ----------------------------------------------------------------------------
# DATASET 1: Blum combined tissue
# ----------------------------------------------------------------------------
# Blum et al. combined tissue dataset — multiple tissue compartments and donors.
# If auto-detection of TCR columns fails, specify overrides:
#   trbv_col  = "<exact_column_name>"
#   clono_col = "<exact_column_name>"
#   sample_col = "<exact_column_name>"
# ----------------------------------------------------------------------------
blum_out <- run_dataset_analysis(
  h5ad_path  = file.path(VAL_DIR, "blum_combined_tissue_data.h5ad"),
  dataset_id = "blum",
  group_vars = "sample_id"
  # trbv_col  = NULL,   # override if needed
  # clono_col = NULL,
  # sample_col = NULL
)

# ----------------------------------------------------------------------------
# DATASET 2: Lechner thyroiditis CD8
# ----------------------------------------------------------------------------
# Lechner et al. — thyroiditis tissue CD8 T cells.
# Only CD8 T cells are expected; a 'subset' column may not exist.
# ----------------------------------------------------------------------------
lechner_out <- run_dataset_analysis(
  h5ad_path  = file.path(VAL_DIR, "lechner_thyroiditis_tissue_cd8_seurat_object.h5ad"),
  dataset_id = "lechner",
  group_vars = "sample_id"
)

# ----------------------------------------------------------------------------
# DATASET 3: Luoma colitis CD8
# ----------------------------------------------------------------------------
# Luoma et al. — colitis tissue CD8 T cells from patients and healthy controls.
# ----------------------------------------------------------------------------
luoma_out <- run_dataset_analysis(
  h5ad_path  = file.path(VAL_DIR, "luoma_colitis_cd8_seurat_object.h5ad"),
  dataset_id = "luoma",
  group_vars = "sample_id"
)


# ── 5. Cross-dataset summary --------------------------------------------------

message("\n", strrep("=", 70))
message("CROSS-DATASET SUMMARY")
message(strrep("=", 70))

# Collect per-dataset correlation results
correlation_summary <- purrr::map_dfr(
  list(blum = blum_out, lechner = lechner_out, luoma = luoma_out),
  function(ds_out) {
    comp <- ds_out$results$comparison
    if (is.null(comp$spearman) || all(is.na(comp$spearman))) {
      return(data.frame(
        spearman_rho_H     = NA_real_, spearman_p_H     = NA_real_,
        pearson_r_H        = NA_real_, pearson_p_H      = NA_real_,
        spearman_rho_clon  = NA_real_, spearman_p_clon  = NA_real_,
        n_groups           = nrow(comp$table)
      ))
    }
    data.frame(
      spearman_rho_H    = round(comp$spearman$estimate,           4),
      spearman_p_H      = round(comp$spearman$p.value,            4),
      pearson_r_H       = round(comp$pearson$estimate,            4),
      pearson_p_H       = round(comp$pearson$p.value,             4),
      spearman_rho_clon = round(comp$clonality_spearman$estimate, 4),
      spearman_p_clon   = round(comp$clonality_spearman$p.value,  4),
      n_groups          = nrow(comp$table)
    )
  },
  .id = "dataset"
)

cat("\n--- H_TRBV vs H_clonotype Spearman correlation per dataset ---\n")
print(correlation_summary, row.names = FALSE)

write.csv(correlation_summary,
          file.path(OUTPUT_DIR, "cross_dataset_correlation_summary.csv"),
          row.names = FALSE)

# Combined diversity tables (one row per sample × dataset)
all_summary <- purrr::map_dfr(
  list(blum = blum_out, lechner = lechner_out, luoma = luoma_out),
  ~ .x$summary_table,
  .id = "dataset"
)

write.csv(all_summary,
          file.path(OUTPUT_DIR, "cross_dataset_diversity_summary.csv"),
          row.names = FALSE)

# ── 6. Cross-dataset scatter: H_TRBV vs H_clonotype --------------------------

# Combine all comparison tables
all_comp_tables <- purrr::map_dfr(
  list(blum = blum_out, lechner = lechner_out, luoma = luoma_out),
  ~ .x$results$comparison$table,
  .id = "dataset"
)

if (nrow(all_comp_tables) >= 3) {

  # Overall correlation across all samples × datasets
  sp_all <- tryCatch(
    cor.test(all_comp_tables$H_TRBV, all_comp_tables$H_clonotype,
             method = "spearman", exact = FALSE),
    error = function(e) NULL
  )
  overall_sub <- if (!is.null(sp_all)) {
    paste0("Overall Spearman rho = ", round(sp_all$estimate, 3),
           "  (p = ", signif(sp_all$p.value, 3), ", n = ", nrow(all_comp_tables), ")")
  } else ""

  p_cross <- ggplot(all_comp_tables,
                    aes(x = H_TRBV, y = H_clonotype,
                        color = dataset, label = sample_id)) +
    geom_point(size = 3, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, color = "grey30",
                linetype = "dashed", linewidth = 0.7) +
    geom_text_repel(size = 2.6, max.overlaps = 15) +
    geom_abline(slope = 1, intercept = 0,
                color = "red", linetype = "dotted", linewidth = 0.6) +
    facet_wrap(vars(dataset), scales = "free") +
    scale_color_manual(
      values = c("blum"    = "#2166AC",
                 "lechner" = "#D6604D",
                 "luoma"   = "#4DAF4A"),
      name = "Dataset"
    ) +
    labs(
      title    = "TRBV Usage Diversity vs True Clonotype Diversity — All Datasets",
      subtitle = overall_sub,
      x        = expression(H[TRBV]~"(TRBV usage diversity proxy)"),
      y        = expression(H[clonotype]~"(true CDR3-level diversity)"),
      caption  = "Red dotted line: 1:1 reference. TRBV proxy cannot resolve CDR3 clonotypes."
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")

  ggsave(file.path(OUTPUT_DIR, "cross_dataset_h_trbv_vs_h_clonotype.pdf"),
         p_cross, width = 13, height = 5)
  message("Cross-dataset scatter saved.")

  # Cross-dataset clonality comparison
  p_cross_clon <- ggplot(all_comp_tables,
                          aes(x = TRBV_clonality_proxy, y = true_clonality,
                              color = dataset, label = sample_id)) +
    geom_point(size = 3, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, color = "grey30",
                linetype = "dashed", linewidth = 0.7) +
    geom_text_repel(size = 2.6, max.overlaps = 15) +
    geom_abline(slope = 1, intercept = 0,
                color = "red", linetype = "dotted", linewidth = 0.6) +
    facet_wrap(vars(dataset), scales = "free") +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(
      values = c("blum" = "#2166AC", "lechner" = "#D6604D", "luoma" = "#4DAF4A"),
      name   = "Dataset"
    ) +
    labs(
      title = "TRBV Clonality Proxy vs True Clonality — All Datasets",
      x     = "TRBV clonality proxy (1 - H_norm_TRBV)",
      y     = "True clonality (1 - H_norm_clonotype)"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")

  ggsave(file.path(OUTPUT_DIR, "cross_dataset_clonality_comparison.pdf"),
         p_cross_clon, width = 13, height = 5)
}


# ── 7. Multi-segment analysis on each dataset --------------------------------
#
# Extends from TRBV-only to TRBJ, VJ-combo, and VDJ-combo proxies, if the
# data contains J-gene information.  Skipped silently if trbj_gene is absent.

message("\n===== Multi-Segment Analysis (where J-gene data is available) =====\n")

run_multiseg_if_available <- function(ds_out, dataset_id) {
  meta <- ds_out$meta
  if (!"trbj_gene" %in% colnames(meta)) {
    message("  [", dataset_id, "] No trbj_gene column — skipping multi-segment analysis.")
    return(invisible(NULL))
  }
  message("  [", dataset_id, "] Running multi-segment analysis...")
  meta_ms <- add_vdj_combos(meta)
  ms <- run_multisegment_analysis(meta_ms, group_vars = "sample_id",
                                   min_cells = MIN_CELLS)
  print_multisegment_summary(ms)

  write.csv(ms$correlations,
            file.path(OUTPUT_DIR, paste0(dataset_id, "_multisegment_correlations.csv")),
            row.names = FALSE)

  p <- plot_correlation_bars(ms$correlations) +
    labs(title = paste("Multi-Segment Correlations —", dataset_id))
  ggsave(file.path(OUTPUT_DIR, paste0(dataset_id, "_multisegment_corr_bars.pdf")),
         p, width = 8, height = 5)

  invisible(ms)
}

run_multiseg_if_available(blum_out,    "blum")
run_multiseg_if_available(lechner_out, "lechner")
run_multiseg_if_available(luoma_out,   "luoma")


# ── 8. Final interpretation notes --------------------------------------------

cat("
=============================================================================
VALIDATION ANALYSIS COMPLETE
=============================================================================

Outputs written to: ", OUTPUT_DIR, "

Key files:
  cross_dataset_correlation_summary.csv   — Spearman rho per dataset
  cross_dataset_diversity_summary.csv     — All samples/groups combined
  cross_dataset_h_trbv_vs_h_clonotype.pdf — Main validation scatter
  <dataset>_tcr_diversity_summary.csv     — Per-sample metrics per dataset

INTERPRETATION:
  A high Spearman rho (> 0.7) between H_TRBV and H_clonotype in a given
  dataset suggests the TRBV usage proxy is a reasonable diversity surrogate
  in that context (e.g., sorted CD8 T cells with little V-gene bias).

  A low rho, or the presence of samples in the 'failure_mode_flags.csv'
  files, indicates contexts where the proxy breaks down and full TCR
  sequencing is essential.

  The TRBV proxy CANNOT distinguish clonal expansion within a TRBV family
  from low V-gene diversity caused by germline/developmental bias.
  Do not interpret H_TRBV as 'true TCR diversity.'
=============================================================================
")
