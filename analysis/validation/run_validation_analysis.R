# =============================================================================
# Validation Analysis — Full TCR Diversity (Beta + Alpha + Paired)
# Published Datasets: Blum Tissue, Lechner, Luoma, Blum PBMC
# =============================================================================
#
# DATASETS:
#   blum_combined_tissue_data.h5ad          — Blum et al. combined tissue T cells
#   lechner_thyroiditis_tissue_cd8_seurat_object.h5ad — Lechner thyroiditis CD8
#   luoma_colitis_cd8_seurat_object.h5ad    — Luoma colitis CD8
#   blum_combined_pbmc_data_with_tcr.h5ad   — Blum et al. PBMC with TCR data
#
# GENE SEGMENTS ANALYSED:
#   Beta chain:  TRBV (50 genes), TRBD (2), TRBJ (13)
#   Alpha chain: TRAV (41 genes), TRAJ (61)
#   Combos:      VJ_beta, VJ_alpha, VA:VB paired, Full 4-segment paired
#
# OUTPUTS (output/validation/):
#   Per dataset:
#     <ds>_tcr_diversity_summary.csv
#     <ds>_full_tcr_correlations.csv
#     <ds>_trbv_stacked_bar.pdf, <ds>_trav_stacked_bar.pdf
#     <ds>_traj_stacked_bar.pdf, <ds>_trbj_stacked_bar.pdf
#     <ds>_clonotype_rank_abundance.pdf
#     <ds>_h_trbv_vs_h_clonotype.pdf, <ds>_clonality_comparison.pdf
#     <ds>_bland_altman.pdf, <ds>_within_trbv_diversity.pdf
#     <ds>_full_tcr_correlation_bars.pdf
#     <ds>_failure_mode_flags.csv
#   Cross-dataset:
#     cross_dataset_correlation_summary.csv
#     cross_dataset_full_tcr_correlations.csv
#     cross_dataset_h_trbv_vs_h_clonotype.pdf
#     cross_dataset_clonality_comparison.pdf
# =============================================================================


# ── 0. Setup ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(stringr)
  library(ggplot2); library(ggrepel); library(scales)
  library(data.table); library(RColorBrewer); library(patchwork); library(pheatmap)
})

REPO_ROOT <- here::here()
source(file.path(REPO_ROOT, "R", "tcr_diversity_functions.R"))

VAL_DIR    <- file.path(REPO_ROOT, "validation")
OUTPUT_DIR <- file.path(REPO_ROOT, "output", "validation")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

MIN_CELLS    <- 30
DOWNSAMPLE_N <- 200

DATASET_COLORS <- c(
  blum_tissue = "#2166AC",
  lechner     = "#D6604D",
  luoma       = "#4DAF4A",
  blum_pbmc   = "#9B59B6"
)


# ── 1. h5ad loading and column normalisation ----------------------------------

load_h5ad_metadata <- function(path) {
  stopifnot(file.exists(path))
  message("  Loading: ", basename(path))
  if (requireNamespace("zellkonverter", quietly = TRUE)) {
    sce  <- zellkonverter::readH5AD(path, X = FALSE, use_hdf5 = FALSE)
    meta <- as.data.frame(SingleCellExperiment::colData(sce))
    return(tibble::rownames_to_column(meta, "barcode_aggr"))
  }
  if (requireNamespace("anndata", quietly = TRUE)) {
    ad   <- anndata::read_h5ad(path)
    meta <- as.data.frame(ad$obs)
    return(tibble::rownames_to_column(meta, "barcode_aggr"))
  }
  stop("Install zellkonverter (BiocManager) or anndata (CRAN).")
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Column detection patterns — beta chain
.TRBV_PAT  <- c("^v_gene$","^trbv$","^trbv_gene$","^v_gene_beta$",
                  "IR_VDJ_1_v_call","^TCR_Beta_V_gene$","^beta_v_gene$","^vb_gene$","v_call")
.TRBJ_PAT  <- c("^j_gene$","^trbj$","^trbj_gene$","^j_gene_beta$",
                  "IR_VDJ_1_j_call","^TCR_Beta_J_gene$","^beta_j_gene$","^jb_gene$","j_call")
.TRBD_PAT  <- c("^d_gene$","^trbd$","^trbd_gene$","^d_gene_beta$",
                  "IR_VDJ_1_d_call","^TCR_Beta_D_gene$","^beta_d_gene$")
# Alpha chain
.TRAV_PAT  <- c("^trav$","^trav_gene$","^v_gene_alpha$",
                  "IR_VJ_1_v_call","^TCR_Alpha_V_gene$","^alpha_v_gene$","^va_gene$")
.TRAJ_PAT  <- c("^traj$","^traj_gene$","^j_gene_alpha$",
                  "IR_VJ_1_j_call","^TCR_Alpha_J_gene$","^alpha_j_gene$","^ja_gene$")
# Clonotype + sample
.CLONO_PAT <- c("^clonotype_id$","^clone_id$","^cloneid$","^cdr3_beta$","^cdr3_b_aa$",
                  "^IR_obs_vj_1_junction_aa$","^TCR_clonotype$","^clonotype$","^clone$",
                  "^raw_clonotype_id$","^ct_id$")
.SAMPLE_PAT <- c("^sample_id$","^sample$","^donor$","^patient$","^patient_id$",
                   "^donor_id$","^orig.ident$","^batch$","^library_id$")

.find_col <- function(nms, pats) {
  for (p in pats) { h <- grep(p, nms, ignore.case=TRUE, value=TRUE); if (length(h)) return(h[1]) }
  NULL
}

#' Auto-detect and normalise all TCR column names
normalise_tcr_metadata <- function(meta, dataset_id,
                                    trbv_col=NULL, trbj_col=NULL, trbd_col=NULL,
                                    trav_col=NULL, traj_col=NULL,
                                    clono_col=NULL, sample_col=NULL) {
  cols <- colnames(meta)

  trbv_col  <- trbv_col  %||% .find_col(cols, .TRBV_PAT)
  trbj_col  <- trbj_col  %||% .find_col(cols, .TRBJ_PAT)
  trbd_col  <- trbd_col  %||% .find_col(cols, .TRBD_PAT)
  trav_col  <- trav_col  %||% .find_col(cols, .TRAV_PAT)
  traj_col  <- traj_col  %||% .find_col(cols, .TRAJ_PAT)
  clono_col <- clono_col %||% .find_col(cols, .CLONO_PAT)
  sample_col<- sample_col%||% .find_col(cols, .SAMPLE_PAT)

  cat(sprintf(
    "\n  [%s] Column mapping:\n    TRBV -> %s | TRBJ -> %s | TRBD -> %s\n    TRAV -> %s | TRAJ -> %s\n    Clonotype -> %s | Sample -> %s\n",
    dataset_id,
    trbv_col %||% "NOT FOUND", trbj_col %||% "NOT FOUND", trbd_col %||% "NOT FOUND",
    trav_col %||% "NOT FOUND", traj_col %||% "NOT FOUND",
    clono_col %||% "NOT FOUND", sample_col %||% paste0("NOT FOUND (using '", dataset_id, "')")
  ))

  rn <- function(df, old, new) {
    if (!is.null(old) && old %in% colnames(df) && old != new) rename(df, !!new := all_of(old))
    else df
  }
  meta <- rn(meta, trbv_col,  "trbv_gene")
  meta <- rn(meta, trbj_col,  "trbj_gene")
  meta <- rn(meta, trbd_col,  "trbd_gene")
  meta <- rn(meta, trav_col,  "trav_gene")
  meta <- rn(meta, traj_col,  "traj_gene")
  meta <- rn(meta, clono_col, "clonotype_id")
  if (!is.null(sample_col) && sample_col != "sample_id")
    meta <- rn(meta, sample_col, "sample_id")
  else if (is.null(sample_col))
    meta$sample_id <- dataset_id

  for (col in c("trbv_gene","trbj_gene","trbd_gene","trav_gene","traj_gene","clonotype_id"))
    if (!col %in% colnames(meta)) meta[[col]] <- NA_character_

  meta$dataset <- dataset_id
  meta
}


# ── 2. Per-dataset analysis wrapper ------------------------------------------

run_dataset_analysis <- function(h5ad_path, dataset_id,
                                  group_vars   = "sample_id",
                                  min_cells    = MIN_CELLS,
                                  downsample_n = DOWNSAMPLE_N,
                                  ...) {

  message("\n", strrep("=", 70))
  message("DATASET: ", dataset_id)
  message(strrep("=", 70))

  meta_raw <- load_h5ad_metadata(h5ad_path)
  cat("Available columns:\n")
  cat(paste(sprintf("  [%03d] %s", seq_along(colnames(meta_raw)), colnames(meta_raw)),
            collapse = "\n"), "\n")

  meta <- normalise_tcr_metadata(meta_raw, dataset_id, ...)

  # Add all chain combo columns
  meta <- add_vdj_combos(meta)           # adds vj_combo (VB:JB) and vdj_combo
  meta <- add_alpha_chain_combos(meta)   # adds vj_alpha, vavb_combo, full_vj_combo

  # Standardise gene calls against canonical lists
  for (cc in list(list("trbv_gene","TRBV"), list("trbj_gene","TRBJ"),
                  list("trbd_gene","TRBD"), list("trav_gene","TRAV"),
                  list("traj_gene","TRAJ"))) {
    if (cc[[1]] %in% colnames(meta))
      meta[[cc[[1]]]] <- standardise_gene_calls(meta[[cc[[1]]]], cc[[2]])
  }
  if ("trav_gene_std" %in% colnames(meta))
    meta$trav_gene_std <- standardise_gene_calls(meta$trav_gene_std, "TRAV")
  if ("traj_gene_std" %in% colnames(meta))
    meta$traj_gene_std <- standardise_gene_calls(meta$traj_gene_std, "TRAJ")

  # Coverage report
  n <- nrow(meta)
  cov <- function(col) {
    if (!col %in% colnames(meta)) return(0L)
    sum(!is.na(meta[[col]]) & meta[[col]] != "" & meta[[col]] != "None")
  }
  cat(sprintf("\n  Total cells:          %d\n", n))
  for (col in c("trbv_gene","trbj_gene","trbd_gene","trav_gene_std","traj_gene_std","clonotype_id")) {
    cat(sprintf("  %-22s: %d (%.1f%%)\n", col, cov(col), 100*cov(col)/n))
  }
  cat(sprintf("  Unique sample_ids:    %d\n", n_distinct(meta$sample_id)))

  # ── Standard TRBV diversity pipeline ─────────────────────────────────────
  results <- run_tcr_diversity_analysis(
    meta, group_vars=group_vars, min_cells=min_cells, downsample_n=downsample_n
  )
  print_correlation_summary(results$comparison)
  summary_tbl <- build_summary_table(results$comparison)
  print(summary_tbl)

  out <- file.path(OUTPUT_DIR, dataset_id)
  write.csv(summary_tbl,
            paste0(out, "_tcr_diversity_summary.csv"), row.names=FALSE)

  bm <- results$benchmark
  if ("failure_mode" %in% colnames(bm$failures))
    write.csv(bm$failures, paste0(out, "_failure_mode_flags.csv"), row.names=FALSE)
  write.csv(bm$within_trbv,
            paste0(out, "_within_trbv_clonotype_details.csv"), row.names=FALSE)

  # ── Full multi-segment analysis (beta + alpha + paired) ───────────────────
  message("\n  Running full multi-segment analysis (beta + alpha + paired)...")
  full_tcr <- tryCatch(
    run_full_tcr_analysis(meta, group_vars=group_vars,
                           min_cells=min_cells, downsample_n=0),
    error = function(e) { warning("Full TCR analysis failed: ", e$message); NULL }
  )
  if (!is.null(full_tcr) && !is.null(full_tcr$correlations))
    write.csv(full_tcr$correlations,
              paste0(out, "_full_tcr_correlations.csv"), row.names=FALSE)

  # ── Plots ─────────────────────────────────────────────────────────────────
  message("  Generating plots...")
  gv1 <- group_vars[1]

  safe_save <- function(expr, path, w=10, h=6) {
    tryCatch({ p <- expr; ggsave(path, p, width=w, height=h) },
             error = function(e) warning(basename(path), " failed: ", e$message))
  }

  safe_save(plot_trbv_stacked_bar(meta, gv1, 20),
            paste0(out, "_trbv_stacked_bar.pdf"))
  safe_save(plot_clonotype_rank_abundance(meta, gv1),
            paste0(out, "_clonotype_rank_abundance.pdf"), w=9)
  safe_save(plot_h_trbv_vs_h_clonotype(results$comparison$table, gv1),
            paste0(out, "_h_trbv_vs_h_clonotype.pdf"), w=7)
  safe_save(plot_clonality_comparison(results$comparison$table, gv1),
            paste0(out, "_clonality_comparison.pdf"), w=7)
  safe_save(plot_bland_altman(results$comparison$table, gv1),
            paste0(out, "_bland_altman.pdf"), w=7)
  safe_save(plot_within_trbv_diversity(bm, gv1),
            paste0(out, "_within_trbv_diversity.pdf"), w=14, h=8)

  # TRBV heatmap
  tryCatch({
    pdf(paste0(out, "_trbv_heatmap.pdf"), width=9, height=8)
    plot_trbv_heatmap(meta, gv1, 25); dev.off()
  }, error = function(e) { dev.off(); warning("trbv_heatmap: ", e$message) })

  # Alpha chain barplots (if data present)
  for (cfg in list(list("trav_gene_std","TRAV"), list("traj_gene_std","TRAJ"),
                   list("trbj_gene","TRBJ"))) {
    col <- cfg[[1]]; nm <- cfg[[2]]
    if (col %in% colnames(meta) &&
        sum(!is.na(meta[[col]]) & meta[[col]] != "") >= min_cells) {
      safe_save(plot_gene_usage_bar(meta, col, gv1, 20,
                                    title=paste(nm, "Gene Usage —", dataset_id)),
                paste0(out, "_", tolower(nm), "_stacked_bar.pdf"))
    }
  }

  # Full TCR correlation bar chart
  if (!is.null(full_tcr) && !is.null(full_tcr$correlations)) {
    safe_save(
      plot_correlation_bars(full_tcr$correlations,
        metric_order=c("TRBV","TRBJ","TRBD","TRAV","TRAJ",
                        "VJ_beta","VJ_alpha","VA_VB","Full_paired")) +
        labs(title=paste("Full TCR Segment Correlations —", dataset_id)),
      paste0(out, "_full_tcr_correlation_bars.pdf"), w=9, h=6
    )
  }

  message("  Outputs saved: ", out, "_*")

  list(meta=meta, results=results, summary_table=summary_tbl,
       benchmark=bm, full_tcr=full_tcr)
}


# ── 3. Run all four datasets --------------------------------------------------

blum_tissue_out <- run_dataset_analysis(
  file.path(VAL_DIR, "blum_combined_tissue_data.h5ad"),
  dataset_id = "blum_tissue",
  group_vars = "sample_id"
)

lechner_out <- run_dataset_analysis(
  file.path(VAL_DIR, "lechner_thyroiditis_tissue_cd8_seurat_object.h5ad"),
  dataset_id = "lechner",
  group_vars = "sample_id"
)

luoma_out <- run_dataset_analysis(
  file.path(VAL_DIR, "luoma_colitis_cd8_seurat_object.h5ad"),
  dataset_id = "luoma",
  group_vars = "sample_id"
)

blum_pbmc_out <- run_dataset_analysis(
  file.path(VAL_DIR, "blum_combined_pbmc_data_with_tcr.h5ad"),
  dataset_id = "blum_pbmc",
  group_vars = "sample_id"
)

all_outputs <- list(
  blum_tissue = blum_tissue_out,
  lechner     = lechner_out,
  luoma       = luoma_out,
  blum_pbmc   = blum_pbmc_out
)


# ── 4. Cross-dataset summaries ------------------------------------------------

message("\n", strrep("=", 70))
message("CROSS-DATASET SUMMARY")
message(strrep("=", 70))

corr_summary <- purrr::map_dfr(all_outputs, function(ds) {
  comp <- ds$results$comparison
  if (is.null(comp$spearman) || all(is.na(comp$spearman)))
    return(tibble(spearman_rho_H=NA_real_, spearman_p_H=NA_real_,
                  pearson_r_H=NA_real_,   spearman_rho_clon=NA_real_,
                  n_groups=nrow(comp$table)))
  tibble(
    spearman_rho_H    = round(comp$spearman$estimate,           3),
    spearman_p_H      = round(comp$spearman$p.value,            4),
    pearson_r_H       = round(comp$pearson$estimate,            3),
    spearman_rho_clon = round(comp$clonality_spearman$estimate, 3),
    n_groups          = nrow(comp$table)
  )
}, .id="dataset")

print(corr_summary)
write.csv(corr_summary,
          file.path(OUTPUT_DIR, "cross_dataset_correlation_summary.csv"),
          row.names=FALSE)

# Cross-dataset full TCR correlation table
full_corr_all <- purrr::map_dfr(all_outputs, function(ds) {
  if (is.null(ds$full_tcr) || is.null(ds$full_tcr$correlations)) return(NULL)
  ds$full_tcr$correlations
}, .id="dataset")

if (nrow(full_corr_all) > 0) {
  write.csv(full_corr_all,
            file.path(OUTPUT_DIR, "cross_dataset_full_tcr_correlations.csv"),
            row.names=FALSE)
  cat("\n--- Full TCR segment correlations (all datasets) ---\n")
  print(full_corr_all %>%
          group_by(metric) %>%
          summarise(mean_rho=round(mean(spearman_rho,na.rm=TRUE),3),
                    sd_rho=round(sd(spearman_rho,na.rm=TRUE),3),
                    n_datasets=sum(!is.na(spearman_rho))) %>%
          arrange(desc(abs(mean_rho))))
}

# Combined diversity table and scatter
all_comp <- purrr::map_dfr(all_outputs, ~ .x$results$comparison$table, .id="dataset")
write.csv(all_comp,
          file.path(OUTPUT_DIR, "cross_dataset_diversity_summary.csv"),
          row.names=FALSE)

if (nrow(all_comp) >= 3) {
  sp_all <- tryCatch(cor.test(all_comp$H_TRBV, all_comp$H_clonotype,
                               method="spearman", exact=FALSE), error=function(e) NULL)
  sub_txt <- if (!is.null(sp_all))
    paste0("Overall Spearman ρ = ", round(sp_all$estimate,3),
           "  (p = ", signif(sp_all$p.value,3), ", n = ", nrow(all_comp), ")")
  else ""

  p_cross <- ggplot(all_comp, aes(H_TRBV, H_clonotype,
                                   color=dataset, label=sample_id)) +
    geom_point(size=3, alpha=.85) +
    geom_smooth(method="lm", se=TRUE, color="grey35",
                linetype="dashed", linewidth=.7) +
    geom_text_repel(size=2.5, max.overlaps=15) +
    geom_abline(slope=1, intercept=0, color="red",
                linetype="dotted", linewidth=.6) +
    facet_wrap(vars(dataset), scales="free") +
    scale_color_manual(values=DATASET_COLORS) +
    labs(title="TRBV Proxy vs True Clonotype Diversity — All Datasets",
         subtitle=sub_txt,
         x=expression(H[TRBV]~"(proxy)"),
         y=expression(H[clonotype]~"(true CDR3 diversity)")) +
    theme_bw(base_size=12) + theme(legend.position="bottom")

  ggsave(file.path(OUTPUT_DIR, "cross_dataset_h_trbv_vs_h_clonotype.pdf"),
         p_cross, width=14, height=6)

  p_clon <- ggplot(all_comp, aes(TRBV_clonality_proxy, true_clonality,
                                   color=dataset, label=sample_id)) +
    geom_point(size=3, alpha=.85) +
    geom_smooth(method="lm", se=TRUE, color="grey35",
                linetype="dashed", linewidth=.7) +
    geom_text_repel(size=2.5, max.overlaps=15) +
    geom_abline(slope=1, intercept=0, color="red",
                linetype="dotted", linewidth=.6) +
    facet_wrap(vars(dataset), scales="free") +
    scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,1)) +
    scale_color_manual(values=DATASET_COLORS) +
    labs(title="Clonality Proxy vs True Clonality — All Datasets",
         x="TRBV clonality proxy", y="True clonality") +
    theme_bw(base_size=12) + theme(legend.position="bottom")

  ggsave(file.path(OUTPUT_DIR, "cross_dataset_clonality_comparison.pdf"),
         p_clon, width=14, height=6)
}

message("\nAll outputs saved to: ", OUTPUT_DIR)
