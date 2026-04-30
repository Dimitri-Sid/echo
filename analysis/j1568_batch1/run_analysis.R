# =============================================================================
# J1568 Batch 1 — TCR Diversity Analysis
# =============================================================================
#
# Proof-of-concept analysis for 8 patient samples with matched:
#   - scVDJ/  per-sample Cell Ranger VDJ outputs
#   - scRNA_aggregated/  Cell Ranger aggr output (34,202 cells × 33,538 genes)
#
# Samples (ordered by gem group in the aggr output, alphabetical):
#   Gem -1: P1207    Gem -2: P1219    Gem -3: P1223    Gem -4: P1226
#   Gem -5: P1227    Gem -6: P1229    Gem -7: P1242    Gem -8: P1243
#
# This script demonstrates three usage modes:
#   A. From raw Cell Ranger outputs (no Seurat needed)
#   B. From a Seurat object (with CD3 expression + annotation filtering)
#   C. From a pre-built standalone metadata data frame
#
# Outputs (saved to output_dir):
#   tcr_diversity_summary.csv       — main diversity metric table
#   trbv_stacked_bar.pdf
#   clonotype_rank_abundance.pdf
#   h_trbv_vs_h_clonotype.pdf
#   clonality_comparison.pdf
#   trbv_heatmap.pdf
#   within_trbv_diversity.pdf
#   bland_altman.pdf
#   within_trbv_clonotype_details.csv
#   failure_mode_flags.csv
# =============================================================================

# ---- 0. Setup ----------------------------------------------------------------

# Install missing packages (run once)
# install.packages(c("dplyr","tidyr","purrr","stringr","ggplot2","ggrepel",
#                    "scales","data.table","Matrix","RColorBrewer","patchwork","pheatmap"))

# Source functions from repo root — works regardless of working directory
REPO_ROOT <- here::here()   # requires 'here' package; or set manually:
# REPO_ROOT <- "/path/to/echo"
source(file.path(REPO_ROOT, "R", "tcr_diversity_functions.R"))

# ── Data paths ────────────────────────────────────────────────────────────────
# Raw Cell Ranger outputs are NOT committed to the repo (see .gitignore).
# Set DATA_DIR to wherever your Cell Ranger outputs live on your machine.
DATA_DIR   <- file.path(REPO_ROOT, "data", "j1568_batch1")   # default location
# DATA_DIR <- "/path/to/your/data"                           # or override here

VDJ_DIR    <- file.path(DATA_DIR, "scVDJ")
AGGR_DIR   <- file.path(DATA_DIR, "scRNA_aggregated")
OUTPUT_DIR <- file.path(REPO_ROOT, "output", "j1568_batch1")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Sample names in gem-group order (must match the order used in cellranger aggr)
SAMPLE_NAMES <- c("P1207", "P1219", "P1223", "P1226",
                   "P1227", "P1229", "P1242", "P1243")

# Minimum cells per sample to include in analysis
MIN_CELLS <- 30

# Optional sensitivity analysis: downsample each sample to this many cells
# (set to 0 to skip)
DOWNSAMPLE_N <- 200


# =============================================================================
# MODE A: FROM RAW CELL RANGER OUTPUTS (no Seurat required)
# =============================================================================

# ---- A1. Load VDJ data -------------------------------------------------------
message("\n===== MODE A: Raw Cell Ranger Outputs =====\n")

# build_tcr_metadata() reads filtered_contig_annotations.csv from each sample,
# keeps productive TRB contigs (one per cell, highest UMI wins),
# corrects gem-group suffixes, and cross-references with aggr barcodes.
meta_raw <- build_tcr_metadata(
  vdj_dir      = VDJ_DIR,
  sample_names = SAMPLE_NAMES,
  aggr_dir     = AGGR_DIR           # NULL to skip aggr barcode cross-reference
)

# Inspect the metadata
glimpse(meta_raw)
cat("Cells per sample:\n")
print(table(meta_raw$sample_id))

# At this stage every cell in meta_raw has a productive TRB contig, which is
# strong evidence for T cell identity. Cells without TCR recovery are excluded.


# ---- A2. Run full diversity analysis -----------------------------------------

# Analysis grouped by sample only
results_by_sample <- run_tcr_diversity_analysis(
  metadata     = meta_raw,
  group_vars   = "sample_id",
  min_cells    = MIN_CELLS,
  downsample_n = DOWNSAMPLE_N
)

# Print correlation summary
print_correlation_summary(results_by_sample$comparison)

# Build and save main summary table
summary_table <- build_summary_table(results_by_sample$comparison)
print(summary_table)
write.csv(summary_table,
          file.path(OUTPUT_DIR, "tcr_diversity_summary.csv"),
          row.names = FALSE)
message("Summary table saved.")


# ---- A3. Examine failure mode benchmarks -------------------------------------

benchmark <- results_by_sample$benchmark

# How many clonotypes share each TRBV gene on average?
cat("\n--- Clonotypes-per-TRBV summary (per sample) ---\n")
print(benchmark$clono_per_trbv)

# Samples with discordant TRBV vs clonotype diversity
if ("failure_mode" %in% colnames(benchmark$failures)) {
  cat("\n--- Failure mode flags ---\n")
  print(benchmark$failures %>%
          select(sample_id, H_TRBV, H_clonotype, discordance, failure_mode) %>%
          arrange(desc(abs(discordance))))
}

# TRBV families with a single dominant clone
cat("\n--- TRBV families with >50% dominated by one clone (top 10 by sample) ---\n")
print(head(benchmark$expanded_within_trbv, 10))

# Save benchmark tables
write.csv(benchmark$within_trbv,
          file.path(OUTPUT_DIR, "within_trbv_clonotype_details.csv"),
          row.names = FALSE)
if ("failure_mode" %in% colnames(benchmark$failures)) {
  write.csv(benchmark$failures,
            file.path(OUTPUT_DIR, "failure_mode_flags.csv"),
            row.names = FALSE)
}


# ---- A4. Downsampled sensitivity analysis ------------------------------------

if (!is.null(results_by_sample$downsampled)) {
  ds_summary <- build_summary_table(results_by_sample$downsampled$comparison)
  cat("\n--- Downsampled (n=", DOWNSAMPLE_N, ") summary ---\n", sep = "")
  print(ds_summary)

  # Compare full vs downsampled Shannon diversity
  compare_ds <- full_join(
    summary_table %>% select(sample_id, H_TRBV, H_clonotype) %>%
      rename(H_TRBV_full = H_TRBV, H_clono_full = H_clonotype),
    ds_summary   %>% select(sample_id, H_TRBV, H_clonotype) %>%
      rename(H_TRBV_ds   = H_TRBV, H_clono_ds   = H_clonotype),
    by = "sample_id"
  )
  cat("\nFull vs downsampled H_clonotype:\n")
  print(compare_ds)
}


# ---- A5. Generate all diagnostic plots ---------------------------------------

message("\nGenerating plots...")

# a. TRBV usage stacked bar
p_bar <- plot_trbv_stacked_bar(meta_raw, group_var = "sample_id", top_n_genes = 20)
ggsave(file.path(OUTPUT_DIR, "trbv_stacked_bar.pdf"),
       p_bar, width = 10, height = 6)

# b. Clonotype rank-abundance curves
p_rank <- plot_clonotype_rank_abundance(meta_raw, group_var = "sample_id")
ggsave(file.path(OUTPUT_DIR, "clonotype_rank_abundance.pdf"),
       p_rank, width = 9, height = 6)

# c. H_TRBV vs H_clonotype scatterplot
p_scatter <- plot_h_trbv_vs_h_clonotype(
  results_by_sample$comparison$table,
  group_var = "sample_id"
)
ggsave(file.path(OUTPUT_DIR, "h_trbv_vs_h_clonotype.pdf"),
       p_scatter, width = 7, height = 6)

# d. Clonality proxy vs true clonality
p_clon <- plot_clonality_comparison(
  results_by_sample$comparison$table,
  group_var = "sample_id"
)
ggsave(file.path(OUTPUT_DIR, "clonality_comparison.pdf"),
       p_clon, width = 7, height = 6)

# e. TRBV gene heatmap
pdf(file.path(OUTPUT_DIR, "trbv_heatmap.pdf"), width = 9, height = 8)
plot_trbv_heatmap(meta_raw, group_var = "sample_id", top_n_genes = 25)
dev.off()

# f. Within-TRBV clonotype diversity (failure mode)
p_within <- plot_within_trbv_diversity(benchmark, group_var = "sample_id")
ggsave(file.path(OUTPUT_DIR, "within_trbv_diversity.pdf"),
       p_within, width = 14, height = 8)

# g. Bland-Altman plot
p_ba <- plot_bland_altman(results_by_sample$comparison$table,
                           group_var = "sample_id")
ggsave(file.path(OUTPUT_DIR, "bland_altman.pdf"),
       p_ba, width = 7, height = 6)

message("All plots saved to: ", OUTPUT_DIR)


# =============================================================================
# MODE A (EXTENDED): MULTI-SEGMENT DIVERSITY — TRBV + TRBJ + TRBD + VJ + VDJ
# =============================================================================
#
# Extends the TRBV-only proxy to include J gene, D gene, and VJ/VDJ gene
# combinations, then ranks each by Spearman rho with true clonotype diversity.
#
# IMPORTANT NOTE ON TRBD IN THIS DATASET:
#   Only ~17% of cells have an assigned D gene (Cell Ranger often cannot
#   unambiguously assign the short TRBD segment). TRBD diversity is reported
#   for completeness but should not be used as a primary proxy. VJ combination
#   diversity is the recommended enhanced proxy when D is unavailable.
# =============================================================================

message("\n===== EXTENDED: Multi-Segment Diversity Analysis =====\n")

# ---- Add VJ/VDJ combination columns to the metadata -------------------------
meta_multi <- add_vdj_combos(meta_raw)

# Quick summary of what is available
cat("Gene-segment data availability:\n")
cat("  TRBV assigned:", sum(!is.na(meta_multi$trbv_gene)), "cells\n")
cat("  TRBJ assigned:", sum(!is.na(meta_multi$trbj_gene)), "cells\n")
cat("  TRBD assigned:", sum(!is.na(meta_multi$trbd_gene)), "cells\n")
cat("  VJ combo:     ", sum(!is.na(meta_multi$vj_combo)),  "cells\n")
cat("  VDJ combo:    ", sum(!is.na(meta_multi$vdj_combo)), "cells\n")
cat("  Unique TRBV:  ", n_distinct(meta_multi$trbv_gene, na.rm = TRUE), "\n")
cat("  Unique TRBJ:  ", n_distinct(meta_multi$trbj_gene, na.rm = TRUE), "\n")
cat("  Unique TRBD:  ", n_distinct(meta_multi$trbd_gene, na.rm = TRUE), "\n")
cat("  Unique VJ:    ", n_distinct(meta_multi$vj_combo,  na.rm = TRUE), "\n")
cat("  Unique VDJ:   ", n_distinct(meta_multi$vdj_combo, na.rm = TRUE), "\n")
cat("  Unique clonotypes:", n_distinct(meta_multi$clonotype_id, na.rm = TRUE), "\n\n")

# ---- Run multi-segment analysis ---------------------------------------------
ms_results <- run_multisegment_analysis(
  metadata     = meta_multi,
  group_vars   = "sample_id",
  min_cells    = MIN_CELLS,
  downsample_n = 0
)

# Print ranked correlation table
print_multisegment_summary(ms_results)

# Save correlation summary
write.csv(ms_results$correlations,
          file.path(OUTPUT_DIR, "multisegment_correlations.csv"),
          row.names = FALSE)

# Save full metrics table
write.csv(ms_results$metrics,
          file.path(OUTPUT_DIR, "multisegment_metrics_wide.csv"),
          row.names = FALSE)

# ---- Plots ------------------------------------------------------------------

# h. Correlation bar chart (ranks all metrics)
p_bars <- plot_correlation_bars(ms_results$correlations)
ggsave(file.path(OUTPUT_DIR, "multisegment_correlation_bars.pdf"),
       p_bars, width = 8, height = 5)

# i. Faceted scatter: each metric vs H_clonotype
p_panel <- plot_metric_correlation_panel(ms_results, group_var = "sample_id")
ggsave(file.path(OUTPUT_DIR, "multisegment_correlation_panel.pdf"),
       p_panel, width = 13, height = 9)

# j. Pairs plot: all H metrics vs each other
p_pairs <- plot_multisegment_pairs(ms_results)
if (!is.null(p_pairs)) {
  ggsave(file.path(OUTPUT_DIR, "multisegment_pairs.pdf"),
         p_pairs, width = 12, height = 10)
}

message("Multi-segment outputs saved to: ", OUTPUT_DIR)


# =============================================================================
# MODE A (COMPARTMENT): CD4 vs CD8 STRATIFIED DIVERSITY
# =============================================================================
#
# The NanoString nCounter TCR Diversity assay is designed for use on sorted or
# defined T cell compartments — not bulk T cells. The correct validation question
# is: within the CD8 compartment, does TRBV (or VJ) diversity track CD8-specific
# clonotype diversity? Ditto for CD4.
#
# This section:
#   1. Classifies all 34,202 cells as CD8 / CD4 / Treg / DN / Non-T using marker
#      gene expression from the aggr matrix (CD8A, CD8B, CD4, CD3D, CD3E, FOXP3).
#   2. Joins annotations onto the VDJ-matched T cells.
#   3. Runs the full diversity analysis stratified by sample_id + subset.
#   4. Compares TRBV proxy vs true clonotype diversity separately for CD4 and CD8.
#   5. Generates all plots faceted by subset.
# =============================================================================

message("\n===== COMPARTMENT ANALYSIS: CD4 vs CD8 Stratification =====\n")

# ---- Annotate T cell subsets from expression matrix -------------------------
message("Annotating T cell subsets from aggr expression matrix...")
subset_annot <- annotate_tcell_subsets_from_aggr(
  aggr_dir        = AGGR_DIR,
  count_threshold = 1,       # >= 1 UMI to call a marker positive
  split_treg      = TRUE     # label Tregs separately from CD4
)

# ---- Join annotations onto VDJ metadata ------------------------------------
message("\nJoining subset annotations onto VDJ metadata...")
meta_subset <- add_subset_annotation(
  vdj_metadata = meta_raw,
  subset_annot = subset_annot,
  keep_non_t   = FALSE   # drop Non-T cells (shouldn't have TCR anyway)
)

cat("\nCells per sample × subset:\n")
print(table(meta_subset$sample_id, meta_subset$subset))

# How many VDJ cells fall into each compartment?
cat("\nProportion of VDJ-matched cells per subset:\n")
print(round(prop.table(table(meta_subset$subset)) * 100, 1))


# ---- Run diversity analysis: sample_id × subset ----------------------------

# Only analyse CD4, CD8, Treg where we have enough cells
meta_for_analysis <- meta_subset %>%
  filter(subset %in% c("CD4", "CD8", "Treg"))

results_by_subset <- run_tcr_diversity_analysis(
  metadata     = meta_for_analysis,
  group_vars   = c("sample_id", "subset"),
  min_cells    = 20,          # lower threshold since subsets have fewer cells
  downsample_n = 0
)

cat("\n--- Compartment diversity summary ---\n")
print(build_summary_table(results_by_subset$comparison))

write.csv(build_summary_table(results_by_subset$comparison),
          file.path(OUTPUT_DIR, "compartment_diversity_summary.csv"),
          row.names = FALSE)


# ---- Correlation within each subset ----------------------------------------

div_tbl <- results_by_subset$comparison$table

# Report correlation per subset separately
for (ss in unique(div_tbl$subset)) {
  sub_tbl <- div_tbl %>% filter(subset == ss)
  cat(sprintf("\n--- %s subset (n = %d samples) ---\n", ss, nrow(sub_tbl)))
  if (nrow(sub_tbl) >= 3) {
    sp <- cor.test(sub_tbl$H_TRBV, sub_tbl$H_clonotype,
                   method = "spearman", exact = FALSE)
    cat(sprintf("  H_TRBV vs H_clonotype: Spearman rho = %.3f (p = %.4f)\n",
                sp$estimate, sp$p.value))
  } else {
    cat("  Insufficient samples for correlation\n")
  }
}


# ---- Run multi-segment analysis by compartment ------------------------------

meta_multi_subset <- add_vdj_combos(meta_for_analysis)

ms_subset <- run_multisegment_analysis(
  metadata   = meta_multi_subset,
  group_vars = c("sample_id", "subset"),
  min_cells  = 20
)

cat("\n--- Multi-segment correlations (all subsets pooled) ---\n")
print_multisegment_summary(ms_subset)

write.csv(ms_subset$correlations,
          file.path(OUTPUT_DIR, "compartment_multisegment_correlations.csv"),
          row.names = FALSE)
write.csv(ms_subset$metrics,
          file.path(OUTPUT_DIR, "compartment_multisegment_metrics.csv"),
          row.names = FALSE)

# Subset-specific correlation table
for (ss in c("CD8", "CD4")) {
  sub_metrics <- ms_subset$metrics %>% filter(subset == ss)
  if (nrow(sub_metrics) < 3) next

  h_cols <- grep("^H_(TRBV|TRBJ|VJ|VDJ)$", colnames(sub_metrics), value = TRUE)
  cat(sprintf("\n--- %s only: rho with H_clonotype ---\n", ss))
  for (hc in h_cols) {
    complete <- !is.na(sub_metrics[[hc]]) & !is.na(sub_metrics$H_clonotype)
    if (sum(complete) < 3) { cat(sprintf("  %-4s: insufficient data\n", sub("H_","",hc))); next }
    sp <- cor.test(sub_metrics[[hc]][complete], sub_metrics$H_clonotype[complete],
                   method = "spearman", exact = FALSE)
    cat(sprintf("  %-4s: rho = %.3f  (p = %.4f, n = %d)\n",
                sub("H_","",hc), sp$estimate, sp$p.value, sum(complete)))
  }
}


# ---- Compartment-stratified plots ------------------------------------------

# k. TRBV stacked bar, faceted by subset
for (ss in c("CD4", "CD8")) {
  sub_meta <- meta_subset %>% filter(subset == ss)
  if (nrow(sub_meta) < 50) next
  p <- plot_trbv_stacked_bar(sub_meta, group_var = "sample_id", top_n_genes = 15) +
    labs(title = paste("TRBV Gene Usage —", ss, "T cells"),
         subtitle = paste0("NanoString-style TRBV usage diversity proxy (", ss, " compartment only)"))
  ggsave(file.path(OUTPUT_DIR, paste0("trbv_stacked_bar_", ss, ".pdf")),
         p, width = 10, height = 6)
}

# l. Multi-segment correlation bars, one per subset
for (ss in c("CD8", "CD4")) {
  sub_corr <- ms_subset$correlations %>% filter(grepl(ss, metric) | !grepl("CD|Treg", metric))

  # Recompute correlations for this subset only
  sub_metrics <- ms_subset$metrics %>% filter(subset == ss)
  if (nrow(sub_metrics) < 3) next

  h_cols <- grep("^H_(TRBV|TRBJ|VJ|VDJ)$", colnames(sub_metrics), value = TRUE)
  sub_cors <- purrr::map_dfr(h_cols, function(hc) {
    label <- sub("^H_", "", hc)
    x <- sub_metrics[[hc]]; y <- sub_metrics$H_clonotype
    ok <- !is.na(x) & !is.na(y)
    if (sum(ok) < 3) return(data.frame(metric=label, n_groups=sum(ok),
                                        spearman_rho=NA_real_, spearman_p=NA_real_))
    sp <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
    data.frame(metric=label, n_groups=sum(ok),
               spearman_rho=round(sp$estimate,3), spearman_p=round(sp$p.value,4))
  })

  p <- plot_correlation_bars(sub_cors) +
    labs(title = paste("Metric Correlation with True Clonotype Diversity —", ss, "T cells"),
         subtitle = paste(ss, "compartment only"))
  ggsave(file.path(OUTPUT_DIR, paste0("multisegment_corr_bars_", ss, ".pdf")),
         p, width = 8, height = 5)
}

# m. Faceted scatter: H_TRBV vs H_clonotype, coloured by subset
p_comp_scatter <- ggplot(
  div_tbl %>% filter(subset %in% c("CD4", "CD8")),
  aes(x = H_TRBV, y = H_clonotype, color = subset, label = sample_id)
) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", linewidth = 0.7) +
  geom_text_repel(size = 2.8, max.overlaps = 12) +
  facet_wrap(vars(subset)) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dotted") +
  scale_color_manual(values = c("CD4" = "#2166AC", "CD8" = "#D6604D")) +
  labs(
    title    = "TRBV Proxy vs True Clonotype Diversity by T Cell Compartment",
    subtitle = "Each point = one sample within one compartment",
    x        = expression(H[TRBV]~"(TRBV usage diversity proxy)"),
    y        = expression(H[clonotype]~"(true CDR3 diversity)")
  ) +
  theme_bw(base_size = 12) + theme(legend.position = "none")

ggsave(file.path(OUTPUT_DIR, "h_trbv_vs_clonotype_by_subset.pdf"),
       p_comp_scatter, width = 10, height = 5)

# n. VJ combo diversity vs H_clonotype, by subset (best proxy)
vj_tbl <- ms_subset$metrics %>%
  filter(subset %in% c("CD4", "CD8"), !is.na(H_VJ), !is.na(H_clonotype))

p_vj_scatter <- ggplot(vj_tbl,
  aes(x = H_VJ, y = H_clonotype, color = subset, label = sample_id)
) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", linewidth = 0.7) +
  geom_text_repel(size = 2.8, max.overlaps = 12) +
  facet_wrap(vars(subset)) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dotted") +
  scale_color_manual(values = c("CD4" = "#2166AC", "CD8" = "#D6604D")) +
  labs(
    title    = "VJ Combo Diversity vs True Clonotype Diversity by Compartment",
    subtitle = "VJ combo is the recommended enhanced proxy (Spearman rho ~0.83 overall)",
    x        = expression(H[VJ]~"(V+J gene combination diversity)"),
    y        = expression(H[clonotype]~"(true CDR3 diversity)")
  ) +
  theme_bw(base_size = 12) + theme(legend.position = "none")

ggsave(file.path(OUTPUT_DIR, "h_vj_vs_clonotype_by_subset.pdf"),
       p_vj_scatter, width = 10, height = 5)

# o. Clonotype rank-abundance by subset
for (ss in c("CD4", "CD8")) {
  sub_meta <- meta_subset %>% filter(subset == ss)
  if (nrow(sub_meta) < 50) next
  p <- plot_clonotype_rank_abundance(sub_meta, group_var = "sample_id") +
    labs(title = paste("Clonotype Rank-Abundance —", ss, "T cells"))
  ggsave(file.path(OUTPUT_DIR, paste0("rank_abundance_", ss, ".pdf")),
         p, width = 9, height = 6)
}

message("\nCompartment analysis outputs saved to: ", OUTPUT_DIR)


# =============================================================================
# MODE B: FROM A SEURAT OBJECT
# =============================================================================
# This section shows how to run the analysis if you have already built a Seurat
# object with cell type annotations (e.g., from Harmony integration + Louvain
# clustering + marker-based annotation).
# =============================================================================

# ---- B0. Build a minimal Seurat object from the aggr matrix ------------------
#
# (Skip this block if you already have a Seurat object.)
#
# This takes a few minutes. For the full pipeline you would typically also run
# NormalizeData, FindVariableFeatures, ScaleData, RunPCA, RunUMAP, and annotate
# clusters with marker genes before proceeding.

if (FALSE) {   # Set to TRUE to run Seurat workflow

  library(Seurat)
  library(Matrix)

  # Read aggr sparse matrix
  mat <- ReadMtx(
    mtx      = file.path(AGGR_DIR, "matrix.mtx"),
    features = file.path(AGGR_DIR, "features.tsv"),
    cells    = file.path(AGGR_DIR, "barcodes.tsv"),
    feature.column = 2   # use gene symbols
  )

  seurat_obj <- CreateSeuratObject(mat, project = "J1568",
                                   min.cells = 3, min.features = 200)

  # Basic QC filtering
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj,
                                                      pattern = "^MT-")
  seurat_obj <- subset(seurat_obj,
                       subset = nFeature_RNA > 200 &
                                nFeature_RNA < 5000 &
                                percent.mt < 20)

  # Standard preprocessing
  seurat_obj <- seurat_obj %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = 3000) %>%
    ScaleData() %>%
    RunPCA(npcs = 30) %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 0.5) %>%
    RunUMAP(dims = 1:20)

  # For demonstration: label all cells as "T cell" (replace with your annotations)
  # In practice you would use marker genes or SingleR/scType for annotation.
  seurat_obj$cell_type <- "T cell"

  # ---- B1. Merge Seurat metadata with VDJ data --------------------------------

  # Load VDJ (must already have been run in Mode A, or run here again)
  # meta_raw <- ...  (see Mode A above)

  meta_seurat <- extract_seurat_metadata(
    seurat_obj       = seurat_obj,
    vdj_data         = meta_raw,
    tcell_annotation = "cell_type",  # filter cells whose annotation contains "T"
    cd3_threshold    = 1.0           # require CD3D+CD3E+CD3G normalised sum >= 1
  )

  # Now run the analysis exactly as in Mode A
  results_seurat <- run_tcr_diversity_analysis(
    metadata     = meta_seurat,
    group_vars   = "sample_id",
    min_cells    = MIN_CELLS,
    downsample_n = 0
  )
  print_correlation_summary(results_seurat$comparison)

  # If cell type annotations are available, stratify by subset
  if ("cell_type" %in% colnames(meta_seurat)) {
    results_by_subset <- run_tcr_diversity_analysis(
      metadata   = meta_seurat %>%
        filter(!is.na(cell_type), cell_type != ""),
      group_vars = c("sample_id", "cell_type"),
      min_cells  = MIN_CELLS
    )

    # Subset-specific clonality plot
    if (nrow(results_by_subset$comparison$table) >= 2) {
      p_subset <- plot_h_trbv_vs_h_clonotype(
        results_by_subset$comparison$table,
        group_var = "cell_type",
        label_col = "sample_id"
      ) + facet_wrap(vars(cell_type))
      ggsave(file.path(OUTPUT_DIR, "h_trbv_vs_clonotype_by_subset.pdf"),
             p_subset, width = 12, height = 8)
    }
  }

}  # END if(FALSE) Seurat block


# =============================================================================
# MODE C: FROM A STANDALONE METADATA DATA FRAME
# =============================================================================
# Use this if you already have a data frame with one row per cell, containing
# at minimum: sample_id, trbv_gene, clonotype_id.
# =============================================================================

# ---- C1. Simulate / load a standalone metadata data frame --------------------
#
# In real use, load your own CSV:
#   meta_standalone <- read.csv("my_cell_metadata.csv")
#
# Required columns:
#   sample_id    — sample identifier
#   trbv_gene    — TRBV gene call (e.g., "TRBV5-1"); NA if not assigned
#   clonotype_id — unique clonotype identifier (e.g., "P1207_clonotype1")
#
# Optional columns (used for stratification):
#   subset       — T cell subset (e.g., "CD4", "CD8", "Treg")
#   condition    — experimental condition / timepoint
#   patient_id   — patient identifier (if samples = donors)

# For this demo, we use the data we already loaded in Mode A:
meta_standalone <- meta_raw %>%
  select(barcode_aggr, sample_id, trbv_gene, clonotype_id,
         cdr3_aa, trav_gene, cdr3_alpha_aa) %>%
  # Add placeholder columns to demonstrate stratified analysis
  mutate(
    subset    = NA_character_,   # replace with actual annotations if available
    condition = NA_character_
  )

# ---- C2. Basic analysis (sample level) --------------------------------------

results_standalone <- run_tcr_diversity_analysis(
  metadata     = meta_standalone,
  group_vars   = "sample_id",
  min_cells    = MIN_CELLS,
  downsample_n = 0
)
print_correlation_summary(results_standalone$comparison)
print(build_summary_table(results_standalone$comparison))


# ---- C3. Subset-stratified analysis ------------------------------------------
# If cell type annotations exist in meta_standalone$subset, run:

if (!all(is.na(meta_standalone$subset))) {
  results_stratified <- run_tcr_diversity_analysis(
    metadata   = meta_standalone %>% filter(!is.na(subset)),
    group_vars = c("sample_id", "subset"),
    min_cells  = MIN_CELLS
  )
  print(build_summary_table(results_stratified$comparison))
}


# ---- C4. Custom per-sample TRBV and clonotype summaries ----------------------

# You can also call the component functions directly for fine-grained control:

trbv_summary <- compute_trbv_diversity(
  metadata  = meta_standalone,
  group_vars = "sample_id",
  min_cells  = MIN_CELLS
)
cat("\n--- TRBV usage diversity (NanoString-style proxy) per sample ---\n")
print(trbv_summary %>%
  select(sample_id, n_cells_trbv, n_trbv_detected,
         H_TRBV, H_norm_TRBV, TRBV_clonality_proxy, top_TRBV_gene))

clono_summary <- compute_clonotype_diversity(
  metadata   = meta_standalone,
  group_vars = "sample_id",
  min_cells  = MIN_CELLS
)
cat("\n--- True clonotype diversity per sample ---\n")
print(clono_summary %>%
  select(sample_id, n_cells_clonotype, n_clonotypes,
         H_clonotype, H_norm_clonotype, true_clonality, top_clone_fraction))

# Direct diversity function calls (single group, vector input):
example_trbv_counts <- table(meta_standalone$trbv_gene[
  meta_standalone$sample_id == "P1207"])

cat("\n--- P1207 TRBV diversity metrics (direct function calls) ---\n")
cat("Shannon H:         ", round(shannon_diversity(as.numeric(example_trbv_counts)), 4), "\n")
cat("Normalized Shannon:", round(normalized_shannon(as.numeric(example_trbv_counts)), 4), "\n")
cat("Simpson diversity: ", round(simpson_diversity(as.numeric(example_trbv_counts)), 4), "\n")
cat("Clonality proxy:   ", round(clonality_from_shannon(as.numeric(example_trbv_counts)), 4), "\n")
cat("Gini coefficient:  ", round(gini_coefficient(as.numeric(example_trbv_counts)), 4), "\n")


# =============================================================================
# INTERPRETATION NOTES
# =============================================================================
cat("
=============================================================================
INTERPRETATION NOTES
=============================================================================

TRBV Usage Diversity (NanoString-style proxy):
  - Measures the breadth and evenness of TCR beta variable gene usage.
  - Analogous to the NanoString nCounter TCR Diversity score.
  - HIGH H_TRBV = many distinct TRBV genes used relatively equally.
  - LOW  H_TRBV = few TRBV genes dominate (but does NOT mean few clonotypes).

True Clonotype Diversity:
  - Measures diversity of unique CDR3 sequences (clonotypes).
  - HIGH H_clonotype = many unique CDR3s used relatively equally (diverse).
  - LOW  H_clonotype = few clonotypes dominate (expanded/clonal repertoire).
  - HIGH true_clonality = highly clonal (oligoclonal expansion).

Why the proxy can disagree with true diversity (failure modes):
  1. HIGH TRBV diversity, LOW clonotype diversity:
     - Multiple TRBV genes are used, but each contains a single dominant clone.
     - Example: vaccine-driven oligoclonal expansion using diverse TRBV genes.
  2. LOW TRBV diversity, HIGH clonotype diversity:
     - Few TRBV genes are used, but each contains many distinct CDR3s.
     - Example: TRBV-biased repertoire (germline or developmental bias) with
       no true clonal expansion.

The TRBV proxy CANNOT distinguish these scenarios. Use full TCR sequencing
(CDR3 clonotypes) for reliable clonality estimation.
=============================================================================
")
