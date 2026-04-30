# ECHO — Evaluator of Clonality Heuristics in Omics

<img width="300" height="300" alt="ChatGPT Image Apr 30, 2026, 04_59_03 PM" src="https://github.com/user-attachments/assets/bbdfdac7-7983-44a1-b7df-9267f1377d3a" />

## Repository Structure

```
echo/
├── R/
│   └── tcr_diversity_functions.R   # All reusable functions (30+ functions)
├── analysis/
│   └── j1568_batch1/
│       └── run_analysis.R          #
├── output/                         # Generated plots and tables (gitignored)
└── README.md
```

---

## Functions Overview

### Core diversity metrics
```r
shannon_diversity(counts)           # Shannon entropy H
normalized_shannon(counts)          # Pielou's evenness H / log(S)
simpson_diversity(counts)           # 1 - D (unbiased)
clonality_from_shannon(counts)      # 1 - normalized Shannon
gini_coefficient(counts)            # Gini inequality index
```

### Per-group diversity computation
```r
compute_trbv_diversity(metadata, group_vars, min_cells)
compute_clonotype_diversity(metadata, group_vars, min_cells)
compute_gene_segment_diversity(metadata, gene_col, group_vars)  # V, D, J, VJ, VDJ
compare_diversity_metrics(trbv_results, clonotype_results)
```

### Multi-segment analysis
```r
add_vdj_combos(metadata)            # adds vj_combo, vdj_combo columns
run_multisegment_analysis(metadata, group_vars)   # TRBV, TRBJ, TRBD, VJ, VDJ
```

### T cell subset annotation (from scRNA-seq)
```r
annotate_tcell_subsets_from_aggr(aggr_dir)   # CD4 / CD8 / Treg / DN from expression
add_subset_annotation(vdj_metadata, subset_annot)
```

### Gene expression program scoring
```r
define_tcell_programs()             # returns 5-program gene lists
score_gene_programs(aggr_dir, programs, barcodes)
aggregate_program_scores(cell_scores, cell_meta, group_vars)
correlate_programs_with_diversity(program_agg, diversity_tbl)
run_program_diversity_analysis(...)  # full pipeline wrapper
```

### Data loading
```r
load_vdj_data(vdj_dir, sample_names)          # Cell Ranger VDJ outputs
load_aggr_barcodes(aggr_dir, sample_names)    # Cell Ranger aggr barcodes
build_tcr_metadata(vdj_dir, sample_names, aggr_dir)
```

### Main wrapper
```r
run_tcr_diversity_analysis(metadata, group_vars, min_cells, downsample_n)
```

### Plots (all return ggplot2 objects or render via pheatmap)
```r
plot_trbv_stacked_bar()
plot_clonotype_rank_abundance()
plot_h_trbv_vs_h_clonotype()
plot_clonality_comparison()
plot_trbv_heatmap()
plot_within_trbv_diversity()
plot_bland_altman()
plot_metric_correlation_panel()
plot_correlation_bars()
plot_multisegment_pairs()
plot_program_correlation_heatmap()
plot_program_scatter_panel()
plot_combined_ranking()
plot_combined_rank_score()
```

---

## Quick Start

### From raw Cell Ranger outputs

```r
source("R/tcr_diversity_functions.R")

# Sample names in the same order as used in cellranger aggr
SAMPLE_NAMES <- c("P1207", "P1219", "P1223", "P1226",
                   "P1227", "P1229", "P1242", "P1243")

# Load VDJ + cross-reference with aggr barcodes
meta <- build_tcr_metadata(
  vdj_dir      = "path/to/scVDJ",
  sample_names = SAMPLE_NAMES,
  aggr_dir     = "path/to/scRNA_aggregated"
)

# Annotate CD4 / CD8 / Treg from expression matrix
subset_annot <- annotate_tcell_subsets_from_aggr("path/to/scRNA_aggregated")
meta         <- add_subset_annotation(meta, subset_annot)

# Run full diversity analysis stratified by compartment
meta_cd4_cd8 <- meta |> filter(subset %in% c("CD4", "CD8"))
results      <- run_tcr_diversity_analysis(
  metadata   = meta_cd4_cd8,
  group_vars = c("sample_id", "subset"),
  min_cells  = 20
)

# Multi-segment analysis (TRBV, TRBJ, VJ, VDJ)
meta_ms <- add_vdj_combos(meta_cd4_cd8)
ms      <- run_multisegment_analysis(meta_ms, c("sample_id", "subset"))

# Gene expression program scoring
prog <- run_program_diversity_analysis(
  aggr_dir      = "path/to/scRNA_aggregated",
  cell_meta     = meta_cd4_cd8,
  diversity_tbl = ms$metrics,
  group_vars    = c("sample_id", "subset")
)
```

### From a Seurat object

```r
# Extract metadata and attach VDJ data
vdj       <- load_vdj_data("path/to/scVDJ", SAMPLE_NAMES)
meta_seu  <- extract_seurat_metadata(
  seurat_obj       = seurat_obj,
  vdj_data         = vdj,
  tcell_annotation = "cell_type",   # column in seurat_obj@meta.data
  cd3_threshold    = 1.0            # require CD3D+CD3E+CD3G >= 1
)
results <- run_tcr_diversity_analysis(meta_seu, "sample_id")
```

### From a standalone metadata data frame

```r
# Minimum required columns: sample_id, trbv_gene, clonotype_id
# Optional: trbj_gene, trbd_gene, subset, condition, patient_id
meta_df <- read.csv("my_metadata.csv")
results <- run_tcr_diversity_analysis(meta_df, "sample_id", min_cells = 30)
```

---

## Dependencies

```r
install.packages(c(
  "dplyr", "tidyr", "purrr", "stringr",
  "ggplot2", "ggrepel", "scales", "patchwork",
  "data.table", "Matrix",
  "RColorBrewer", "pheatmap"
))
```

---

## Citation

If you use this framework, please cite: *(manuscript in preparation)*
