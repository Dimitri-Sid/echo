# =============================================================================
# Extract real analysis outputs for the validation website (rhdf5 reader)
# Writes CSVs to analysis/validation/data/ — loaded directly by index.html
# =============================================================================
# Dataset column map:
#   Lechner:    TRBV=v_gene_TRB  TRAV=v_gene_TRA  clono=raw_clonotype_id  sample=orig.ident   subset=var.CellType
#   Luoma:      TRBV=v_gene_TRB  TRAV=v_gene_TRA  clono=raw_clonotype_id  sample=sample_id    subset=(CD8 only)
#   Blum PBMC:  TRBV=TRB_v_gene  TRAV=TRA_v_gene  clono=TRB_cdr3          sample=sample_id    subset=cluster_name
#   Blum Tissue:no TCR columns — RNA-seq only; lineage_names + donor
# =============================================================================

suppressPackageStartupMessages({
  library(rhdf5); library(dplyr); library(tidyr); library(purrr)
  library(stringr); library(data.table); library(here)
})

REPO_ROOT <- here::here()
source(file.path(REPO_ROOT, "R", "tcr_diversity_functions.R"))

VAL_DIR  <- file.path(REPO_ROOT, "validation")
DATA_DIR <- file.path(REPO_ROOT, "analysis", "validation", "data")
dir.create(DATA_DIR, showWarnings=FALSE, recursive=TRUE)

MIN_CELLS <- 10


# ── rhdf5 obs reader ----------------------------------------------------------
# h5ad categoricals are stored as integer codes + a __categories lookup.
# This function returns a plain data.frame.

read_h5ad_obs <- function(path) {
  message("  Reading obs from: ", basename(path))
  info <- h5ls(path, recursive=2)
  obs_nodes  <- info[info$group == "/obs", ]
  cats_avail <- info$name[info$group == "/obs/__categories"]

  # All obs column names (both datasets and groups, excluding specials)
  all_obs_cols <- obs_nodes$name[!obs_nodes$name %in% c("__categories","_index","barcodekey")]

  decode_raw <- function(raw, col) {
    # Format A: named list with 'categories' + 'codes' (new pandas h5ad format)
    if (is.list(raw) && all(c("categories","codes") %in% names(raw))) {
      cats  <- as.character(raw$categories)
      codes <- as.integer(raw$codes)
      result <- ifelse(codes < 0L, NA_character_, cats[codes + 1L])
      return(result)
    }
    # Format B: unnamed list of length 2 (older format)
    if (is.list(raw) && length(raw) == 2) {
      cats  <- as.character(raw[[1]])
      codes <- as.integer(raw[[2]])
      result <- ifelse(codes < 0L, NA_character_, cats[codes + 1L])
      return(result)
    }
    # Format C: integer codes + /obs/__categories/<col> lookup
    if (col %in% cats_avail && (is.integer(raw) || is.numeric(raw))) {
      cats <- tryCatch(h5read(path, paste0("/obs/__categories/", col)), error=function(e) NULL)
      if (!is.null(cats)) {
        codes  <- as.integer(raw)
        return(ifelse(codes < 0L, NA_character_, as.character(cats[codes + 1L])))
      }
    }
    # Format D: plain string or numeric
    as.character(raw)
  }

  obs <- list()
  for (col in all_obs_cols) {
    raw <- tryCatch(h5read(path, paste0("/obs/", col)), error=function(e) NULL)
    if (is.null(raw)) next
    decoded <- tryCatch(decode_raw(raw, col), error=function(e) NULL)
    if (!is.null(decoded)) obs[[col]] <- decoded
  }

  # Keep only columns with consistent length
  n <- if (length(obs) > 0) max(lengths(obs)) else 0L
  obs <- obs[lengths(obs) == n]

  # Barcode
  idx <- tryCatch(as.character(h5read(path, "/obs/_index")),     error=function(e) NULL)
  if (is.null(idx)) idx <- tryCatch(as.character(h5read(path, "/obs/barcodekey")), error=function(e) NULL)
  if (is.null(idx)) idx <- as.character(seq_len(n))
  obs[["barcode_aggr"]] <- idx

  as.data.frame(obs, stringsAsFactors=FALSE)
}


# ── Dataset definitions -------------------------------------------------------

# Only Lechner (2,782 cells — smallest dataset)
DATASETS <- list(
  lechner = list(
    file       = "lechner_thyroiditis_tissue_cd8_seurat_object.h5ad",
    trbv_col   = "v_gene_TRB",
    trbj_col   = "j_gene_TRB",
    trbd_col   = "d_gene_TRB",
    trav_col   = "v_gene_TRA",
    traj_col   = "j_gene_TRA",
    clono_col  = "raw_clonotype_id",
    sample_col = "orig.ident",
    subset_col = NULL,           # all CD8 — cluster numbers not cell types
    fixed_subset = "CD8",
    label      = "Lechner (Thyroiditis CD8)"
  )
)


# ── Load and prep one dataset -------------------------------------------------

load_and_prep <- function(ds_id, ds) {
  path <- file.path(VAL_DIR, ds$file)
  if (!file.exists(path)) { message("  File not found: ", path); return(NULL) }

  obs <- tryCatch(read_h5ad_obs(path), error=function(e) {
    message("  Error reading: ", e$message); return(NULL)
  })
  if (is.null(obs)) return(NULL)
  message("  Loaded ", nrow(obs), " cells, ", ncol(obs), " columns")

  # Rename to standard schema
  rn <- function(df, old, new) {
    if (!is.null(old) && old %in% colnames(df) && old != new)
      dplyr::rename(df, !!new := dplyr::all_of(old))
    else df
  }
  obs <- rn(obs, ds$trbv_col,   "trbv_gene")
  obs <- rn(obs, ds$trbj_col,   "trbj_gene")
  obs <- rn(obs, ds$trbd_col,   "trbd_gene")
  obs <- rn(obs, ds$trav_col,   "trav_gene")
  obs <- rn(obs, ds$traj_col,   "traj_gene")
  obs <- rn(obs, ds$clono_col,  "clonotype_id")
  obs <- rn(obs, ds$sample_col, "sample_id")

  for (col in c("trbv_gene","trbj_gene","trbd_gene","trav_gene","traj_gene","clonotype_id"))
    if (!col %in% colnames(obs)) obs[[col]] <- NA_character_

  # Standardise gene calls
  for (cc in list(list("trbv_gene","TRBV"), list("trbj_gene","TRBJ"),
                  list("trbd_gene","TRBD"), list("trav_gene","TRAV"),
                  list("traj_gene","TRAJ")))
    if (!is.null(ds[[paste0(sub("_gene","_col",cc[[1]]))]])  &&
        cc[[1]] %in% colnames(obs))
      obs[[cc[[1]]]] <- standardise_gene_calls(obs[[cc[[1]]]], cc[[2]])

  # Subset annotation
  if (!is.null(ds$fixed_subset)) {
    obs$subset <- ds$fixed_subset
    message("  Fixed subset: all cells = ", ds$fixed_subset)
  } else if (!is.null(ds$subset_col) && ds$subset_col %in% colnames(obs)) {
    obs$subset <- standardise_subset_extended(obs[[ds$subset_col]])
    obs <- obs[obs$subset != "Non_T" & obs$subset != "Unknown", ]
    message("  After subset filter: ", nrow(obs), " cells")
    print(table(obs$subset))
  } else {
    obs$subset <- "Unknown"
  }

  # Add combo columns
  obs <- add_vdj_combos(obs)
  obs <- add_alpha_chain_combos(obs)
  if ("trav_gene_std" %in% colnames(obs))
    obs$trav_gene_std <- standardise_gene_calls(obs$trav_gene_std, "TRAV")
  if ("traj_gene_std" %in% colnames(obs))
    obs$traj_gene_std <- standardise_gene_calls(obs$traj_gene_std, "TRAJ")

  obs$dataset <- ds_id
  obs
}


# ── Extract all outputs -------------------------------------------------------

all_cov   <- list(); all_trbv  <- list(); all_trav  <- list(); all_traj  <- list()
all_div_w <- list(); all_cor_w <- list(); all_div_s <- list(); all_cor_s <- list()
all_trbv_s<- list(); all_segs  <- list()
all_pseudo_pl <- list(); all_pseudo_sm <- list()
all_ra    <- list()

for (ds_id in names(DATASETS)) {
  ds   <- DATASETS[[ds_id]]
  message("\n", strrep("=",60), "\n", ds_id, "\n", strrep("=",60))
  obs  <- load_and_prep(ds_id, ds)
  if (is.null(obs)) next

  has_tcr <- !is.null(ds$trbv_col) &&
             sum(!is.na(obs$trbv_gene) & obs$trbv_gene != "") > 0
  has_clono <- sum(!is.na(obs$clonotype_id) & obs$clonotype_id != "") > 0
  has_sub   <- "subset" %in% colnames(obs) &&
               any(c("CD4","CD8","DP","DN","Treg") %in% unique(obs$subset))
  n <- nrow(obs)

  # Coverage
  all_cov[[ds_id]] <- tibble(
    dataset=ds_id, label=ds$label, n_cells=n,
    n_samples=n_distinct(obs$sample_id),
    n_trbv  = sum(!is.na(obs$trbv_gene) & obs$trbv_gene != ""),
    n_trav  = if ("trav_gene_std" %in% colnames(obs))
                sum(!is.na(obs$trav_gene_std) & obs$trav_gene_std != "") else 0L,
    n_clono = sum(!is.na(obs$clonotype_id) & obs$clonotype_id != ""),
    n_unique_clonotypes = n_distinct(obs$clonotype_id[!is.na(obs$clonotype_id) &
                                                        obs$clonotype_id != ""])
  )

  # Gene usage
  for (cfg in list(
    list(col="trbv_gene",     key="trbv",  out_col="trbv_gene"),
    list(col="trav_gene_std", key="trav",  out_col="trav_gene"),
    list(col="traj_gene_std", key="traj",  out_col="traj_gene")
  )) {
    if (!cfg$col %in% colnames(obs)) next
    valid <- obs %>% filter(!is.na(.data[[cfg$col]]), .data[[cfg$col]] != "")
    if (nrow(valid) < 5) next
    usage <- valid %>%
      count(sample_id, gene=.data[[cfg$col]]) %>%
      group_by(sample_id) %>% mutate(freq=n/sum(n)) %>% ungroup() %>%
      mutate(dataset=ds_id, chain=toupper(cfg$key))
    if (cfg$key == "trbv") all_trbv[[ds_id]] <- usage
    if (cfg$key == "trav") all_trav[[ds_id]] <- usage
    if (cfg$key == "traj") all_traj[[ds_id]] <- usage
  }

  if (!has_tcr || !has_clono) {
    message("  Skipping diversity analysis (no TCR/clonotype data)")
    next
  }

  # Whole-sample diversity
  res_w <- tryCatch(
    run_tcr_diversity_analysis(obs, "sample_id", MIN_CELLS, downsample_n=0),
    error=function(e) { message("  Whole diversity: ",e$message); NULL }
  )
  if (!is.null(res_w)) {
    sum_w <- build_summary_table(res_w$comparison) %>% mutate(dataset=ds_id)
    all_div_w[[ds_id]] <- sum_w
    comp <- res_w$comparison
    sp   <- comp$spearman; pe <- comp$pearson; sc <- comp$clonality_spearman
    is_valid <- function(x) is.list(x) && !is.null(x$estimate)
    all_cor_w[[ds_id]] <- tibble(
      dataset=ds_id, subset="All",
      spearman_rho = if (is_valid(sp)) round(sp$estimate,3) else NA_real_,
      spearman_p   = if (is_valid(sp)) round(sp$p.value,  4) else NA_real_,
      pearson_r    = if (is_valid(pe)) round(pe$estimate, 3) else NA_real_,
      rho_clonality= if (is_valid(sc)) round(sc$estimate, 3) else NA_real_,
      n_groups     = nrow(comp$table)
    )
  }

  # Subset-stratified diversity
  if (has_sub) {
    obs_s <- obs %>% filter(subset %in% c("CD4","CD8","DP","DN","Treg"))
    res_s <- tryCatch(
      run_tcr_diversity_analysis(obs_s, c("sample_id","subset"), MIN_CELLS, downsample_n=0),
      error=function(e) { message("  Subset diversity: ",e$message); NULL }
    )
    if (!is.null(res_s)) {
      all_div_s[[ds_id]] <- build_summary_table(res_s$comparison) %>% mutate(dataset=ds_id)
      div_tbl <- res_s$comparison$table
      all_cor_s[[ds_id]] <- purrr::map_dfr(
        intersect(c("CD4","CD8","DP","DN","Treg"), unique(div_tbl$subset)), function(ss) {
          st <- div_tbl %>% filter(subset==ss)
          if (nrow(st) < 3) return(NULL)
          sp2 <- tryCatch(cor.test(st$H_TRBV, st$H_clonotype, method="spearman",exact=FALSE),
                          error=function(e) NULL)
          if (is.null(sp2)) return(NULL)
          tibble(dataset=ds_id, subset=ss,
                 spearman_rho=round(sp2$estimate,3),
                 spearman_p=round(sp2$p.value,4), n_groups=nrow(st))
        }
      )
      # TRBV usage per subset
      all_trbv_s[[ds_id]] <- obs_s %>%
        filter(!is.na(trbv_gene), trbv_gene != "") %>%
        count(sample_id, subset, gene=trbv_gene) %>%
        group_by(sample_id, subset) %>% mutate(freq=n/sum(n)) %>% ungroup() %>%
        mutate(dataset=ds_id)
    }
  }

  # Full multi-segment correlations
  seg <- tryCatch(
    run_full_tcr_analysis(obs, "sample_id", MIN_CELLS),
    error=function(e) { message("  Full TCR: ",e$message); NULL }
  )
  if (!is.null(seg) && !is.null(seg$correlations))
    all_segs[[ds_id]] <- seg$correlations %>% mutate(dataset=ds_id)

  # Pseudoclonotype analysis
  ps_gv  <- if (has_sub) c("sample_id","subset") else "sample_id"
  ps_obs <- if (has_sub) obs %>% filter(subset %in% c("CD4","CD8","DP","DN","Treg")) else obs
  pseudo <- tryCatch(
    run_pseudoclonotype_analysis(ps_obs, ps_gv, MIN_CELLS),
    error=function(e) { message("  Pseudo: ",e$message); NULL }
  )
  if (!is.null(pseudo) && nrow(pseudo$per_level) > 0) {
    all_pseudo_pl[[ds_id]] <- pseudo$per_level %>% mutate(dataset=ds_id)
    all_pseudo_sm[[ds_id]] <- pseudo$summary   %>% mutate(dataset=ds_id)
  }

  # Rank abundance (top 100 clonotypes per sample)
  all_ra[[ds_id]] <- obs %>%
    filter(!is.na(clonotype_id), clonotype_id != "") %>%
    group_by(sample_id, clonotype_id) %>% summarise(n=n(),.groups="drop") %>%
    group_by(sample_id) %>% arrange(desc(n)) %>%
    mutate(rank=row_number(), freq=n/sum(n)) %>%
    filter(rank <= 100) %>% ungroup() %>%
    mutate(dataset=ds_id)
}


# ── Write CSVs ----------------------------------------------------------------

write_csv_safe <- function(lst, filename) {
  rows <- purrr::compact(lst)
  if (length(rows) == 0) { message("No data: ", filename); return(invisible(NULL)) }
  df <- bind_rows(rows)
  write.csv(df, file.path(DATA_DIR, filename), row.names=FALSE)
  message("Wrote ", nrow(df), " rows -> ", filename)
}

write_csv_safe(all_cov,       "coverage.csv")
write_csv_safe(all_trbv,      "trbv_usage.csv")
write_csv_safe(all_trav,      "trav_usage.csv")
write_csv_safe(all_traj,      "traj_usage.csv")
write_csv_safe(all_div_w,     "diversity_whole.csv")
write_csv_safe(all_cor_w,     "correlation_whole.csv")
write_csv_safe(all_div_s,     "diversity_subset.csv")
write_csv_safe(all_cor_s,     "correlation_subset.csv")
write_csv_safe(all_trbv_s,    "trbv_usage_subset.csv")
write_csv_safe(all_segs,      "segment_correlations.csv")
write_csv_safe(all_pseudo_pl, "pseudo_per_level.csv")
write_csv_safe(all_pseudo_sm, "pseudo_summary.csv")
write_csv_safe(all_ra,        "rank_abundance.csv")

message("\nDone. CSVs in: ", DATA_DIR)
