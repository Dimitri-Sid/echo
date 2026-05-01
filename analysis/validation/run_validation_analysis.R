# =============================================================================
# Validation Analysis — Full TCR Diversity
# RNA-seq vs TCR-seq Gene Usage · CD4/CD8/DP/DN Stratification
# Pseudoclonotype Recapitulation · Datasets: Blum Tissue, Lechner, Luoma, Blum PBMC
# =============================================================================

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
  blum_tissue="#2166AC", lechner="#D6604D",
  luoma="#4DAF4A",       blum_pbmc="#9B59B6"
)
SUBSET_COLORS <- c(
  CD8="#D6604D", CD4="#2166AC", Treg="#4DAF4A",
  DP="#FF7F00",  DN="#9B59B6",  Unknown="grey60"
)
ALL_SUBSETS <- c("CD8","CD4","DP","DN","Treg")


# ── 1. Loading helpers --------------------------------------------------------

`%||%` <- function(a, b) if (!is.null(a)) a else b

.TRBV_PAT   <- c("^v_gene$","^trbv$","^trbv_gene$","^v_gene_beta$","IR_VDJ_1_v_call","v_call")
.TRBJ_PAT   <- c("^j_gene$","^trbj$","^trbj_gene$","IR_VDJ_1_j_call","j_call")
.TRBD_PAT   <- c("^d_gene$","^trbd$","^trbd_gene$","IR_VDJ_1_d_call")
.TRAV_PAT   <- c("^trav$","^trav_gene$","^v_gene_alpha$","IR_VJ_1_v_call","^alpha_v_gene$")
.TRAJ_PAT   <- c("^traj$","^traj_gene$","^j_gene_alpha$","IR_VJ_1_j_call")
.CLONO_PAT  <- c("^clonotype_id$","^clone_id$","^cdr3_beta$","^IR_obs_vj_1_junction_aa$",
                   "^clonotype$","^clone$","^raw_clonotype_id$")
.SAMPLE_PAT <- c("^sample_id$","^sample$","^donor$","^patient$","^patient_id$",
                   "^donor_id$","^orig.ident$","^batch$")

.find_col <- function(nms, pats) {
  for (p in pats) { h <- grep(p, nms, ignore.case=TRUE, value=TRUE); if (length(h)) return(h[1]) }
  NULL
}

normalise_tcr_metadata <- function(meta, dataset_id, ...) {
  overrides <- list(...); cols <- colnames(meta)
  get <- function(key, pats) overrides[[key]] %||% .find_col(cols, pats)

  trbv <- get("trbv_col",.TRBV_PAT); trbj <- get("trbj_col",.TRBJ_PAT)
  trbd <- get("trbd_col",.TRBD_PAT); trav <- get("trav_col",.TRAV_PAT)
  traj <- get("traj_col",.TRAJ_PAT); clono<- get("clono_col",.CLONO_PAT)
  samp <- get("sample_col",.SAMPLE_PAT)

  cat(sprintf(
    "\n  [%s]\n    TRBV->%s TRBJ->%s TRBD->%s\n    TRAV->%s TRAJ->%s\n    Clono->%s Sample->%s\n",
    dataset_id, trbv%||%"?", trbj%||%"?", trbd%||%"?",
    trav%||%"?", traj%||%"?", clono%||%"?", samp%||%paste0("('",dataset_id,"')")
  ))

  rn <- function(df, old, new) {
    if (!is.null(old) && old %in% colnames(df) && old != new)
      rename(df, !!new := all_of(old)) else df
  }
  meta <- rn(meta,trbv,"trbv_gene"); meta <- rn(meta,trbj,"trbj_gene")
  meta <- rn(meta,trbd,"trbd_gene"); meta <- rn(meta,trav,"trav_gene")
  meta <- rn(meta,traj,"traj_gene"); meta <- rn(meta,clono,"clonotype_id")
  if (!is.null(samp) && samp != "sample_id") meta <- rn(meta,samp,"sample_id")
  else if (is.null(samp)) meta$sample_id <- dataset_id

  for (col in c("trbv_gene","trbj_gene","trbd_gene","trav_gene","traj_gene","clonotype_id"))
    if (!col %in% colnames(meta)) meta[[col]] <- NA_character_

  meta$dataset <- dataset_id
  meta
}


# ── 2. Per-dataset analysis wrapper ------------------------------------------

run_dataset_analysis <- function(h5ad_path, dataset_id,
                                  min_cells    = MIN_CELLS,
                                  downsample_n = DOWNSAMPLE_N, ...) {

  message("\n", strrep("=",70), "\nDATASET: ", dataset_id, "\n", strrep("=",70))
  out <- file.path(OUTPUT_DIR, dataset_id)

  # ── Load metadata only first (fast) ─────────────────────────────────────────
  tcr_genes_all <- unlist(TCR_GENES, use.names=FALSE)

  # Load WITH expression matrix (for RNA-seq gene usage)
  message("\n  Loading h5ad with expression matrix...")
  expr_data <- tryCatch(
    load_h5ad_with_expression(h5ad_path, genes_keep=tcr_genes_all),
    error = function(e) {
      warning("Could not load expression: ", e$message,
              "\nFalling back to metadata-only.")
      # Fallback: metadata-only load
      if (requireNamespace("zellkonverter", quietly=TRUE)) {
        sce  <- zellkonverter::readH5AD(h5ad_path, X=FALSE, use_hdf5=FALSE)
        meta <- tibble::rownames_to_column(
          as.data.frame(SingleCellExperiment::colData(sce)), "barcode_aggr")
        list(meta=meta, expr=NULL, genes=character(0))
      } else {
        ad   <- anndata::read_h5ad(h5ad_path)
        meta <- tibble::rownames_to_column(as.data.frame(ad$obs), "barcode_aggr")
        list(meta=meta, expr=NULL, genes=character(0))
      }
    }
  )

  meta <- normalise_tcr_metadata(expr_data$meta, dataset_id, ...)

  # Add all chain combo columns and standardise gene calls
  meta <- add_vdj_combos(meta)
  meta <- add_alpha_chain_combos(meta)
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

  # ── RNA-seq gene calls ───────────────────────────────────────────────────────
  has_expr <- !is.null(expr_data$expr) && length(expr_data$genes) > 0
  if (has_expr) {
    message("\n  Assigning RNA-seq dominant TCR gene calls...")
    meta <- add_rnaseq_gene_calls(meta, expr_data, min_expr=0)
  } else {
    for (col in c("trbv_rna","trav_rna","trbj_rna","traj_rna"))
      meta[[col]] <- NA_character_
    message("  No expression matrix — RNA-seq gene calls skipped.")
  }

  # ── Extended CD4/CD8/DP/DN subset annotation ──────────────────────────────
  message("\n  Detecting T cell subsets (CD4/CD8/DP/DN/Treg)...")
  meta <- add_subset_extended(meta, keep_nont=FALSE, keep_unknown=FALSE)
  has_subset <- any(ALL_SUBSETS %in% unique(meta$subset))

  cat("  Cells per subset:\n")
  print(table(meta$subset))
  cat(sprintf("  Cells total (post-filter): %d\n\n", nrow(meta)))

  safe <- function(expr_plot, path, w=10, h=6) {
    tryCatch({ ggsave(path, expr_plot, width=w, height=h) },
             error=function(e) warning(basename(path),": ",e$message))
  }

  # ── [A] Whole-sample diversity ────────────────────────────────────────────
  message("  [A] Whole-sample diversity analysis...")
  res_whole <- run_tcr_diversity_analysis(
    meta, group_vars="sample_id", min_cells=min_cells, downsample_n=downsample_n
  )
  print_correlation_summary(res_whole$comparison)
  sum_whole <- build_summary_table(res_whole$comparison)
  write.csv(sum_whole, paste0(out,"_tcr_diversity_summary.csv"), row.names=FALSE)

  bm <- res_whole$benchmark
  if ("failure_mode" %in% colnames(bm$failures))
    write.csv(bm$failures, paste0(out,"_failure_mode_flags.csv"), row.names=FALSE)

  # ── [B] CD4/CD8/DP/DN stratified diversity ────────────────────────────────
  res_subset <- NULL; sum_subset <- NULL
  if (has_subset) {
    message("\n  [B] Stratified diversity (sample × subset)...")
    meta_sub <- meta %>% filter(subset %in% ALL_SUBSETS)
    res_subset <- tryCatch(
      run_tcr_diversity_analysis(
        meta_sub, group_vars=c("sample_id","subset"),
        min_cells=20, downsample_n=0
      ),
      error=function(e) { warning("Subset analysis: ",e$message); NULL }
    )
    if (!is.null(res_subset)) {
      sum_subset <- build_summary_table(res_subset$comparison)
      write.csv(sum_subset, paste0(out,"_subset_diversity_summary.csv"), row.names=FALSE)

      div_tbl <- res_subset$comparison$table
      for (ss in intersect(ALL_SUBSETS, unique(div_tbl$subset))) {
        stbl <- div_tbl %>% filter(subset==ss)
        if (nrow(stbl) >= 3) {
          sp <- cor.test(stbl$H_TRBV, stbl$H_clonotype, method="spearman", exact=FALSE)
          message(sprintf("  %-6s: rho=%.3f p=%.4f n=%d", ss, sp$estimate, sp$p.value, nrow(stbl)))
        }
      }

      # Stratified scatter
      if (nrow(div_tbl) >= 2) {
        p_sub <- ggplot(div_tbl, aes(H_TRBV, H_clonotype, color=subset, label=sample_id)) +
          geom_point(size=3) +
          geom_smooth(method="lm",se=TRUE,linetype="dashed",linewidth=.7) +
          geom_text_repel(size=2.5, max.overlaps=12) +
          geom_abline(slope=1, intercept=0, color="red", linetype="dotted") +
          facet_wrap(vars(subset), scales="free") +
          scale_color_manual(values=SUBSET_COLORS) +
          labs(title=paste("TRBV Proxy vs True Diversity by Subset —",dataset_id),
               x=expression(H[TRBV]), y=expression(H[clonotype])) +
          theme_bw(base_size=11) + theme(legend.position="none")
        safe(p_sub, paste0(out,"_subset_scatter.pdf"), w=13, h=6)
      }

      # Rank-abundance per subset
      for (ss in intersect(ALL_SUBSETS, unique(meta_sub$subset))) {
        sub_m <- meta_sub %>% filter(subset==ss)
        if (nrow(sub_m) < min_cells) next
        p_ra <- plot_clonotype_rank_abundance(sub_m,"sample_id") +
          labs(title=paste("Rank-Abundance —",ss,"T cells —",dataset_id))
        safe(p_ra, paste0(out,"_rank_abundance_",ss,".pdf"), w=9)
      }
    }
  }

  # ── [C] Full multi-segment analysis ─────────────────────────────────────
  message("\n  [C] Full multi-segment TCR analysis (beta + alpha + paired)...")
  full_tcr <- tryCatch(
    run_full_tcr_analysis(meta, group_vars="sample_id", min_cells=min_cells),
    error=function(e) { warning(e$message); NULL }
  )
  if (!is.null(full_tcr$correlations))
    write.csv(full_tcr$correlations, paste0(out,"_full_tcr_correlations.csv"), row.names=FALSE)

  # ── [D] Pseudoclonotype analysis (stratified by subset) ──────────────────
  message("\n  [D] Pseudoclonotype recapitulation analysis...")
  pseudo_gv <- if (has_subset) c("sample_id","subset") else "sample_id"
  meta_ps   <- if (has_subset) meta %>% filter(subset %in% ALL_SUBSETS) else meta

  pseudo <- tryCatch(
    run_pseudoclonotype_analysis(meta_ps, group_vars=pseudo_gv, min_cells=min_cells),
    error=function(e) { warning(e$message); NULL }
  )
  if (!is.null(pseudo) && nrow(pseudo$per_level) > 0) {
    write.csv(pseudo$per_level, paste0(out,"_pseudoclonotype_metrics.csv"), row.names=FALSE)
    write.csv(pseudo$summary,   paste0(out,"_pseudoclonotype_summary.csv"), row.names=FALSE)
  }

  # ── [E] RNA-seq vs TCR-seq comparison ────────────────────────────────────
  rna_comparison <- list()
  if (has_expr && "trbv_rna" %in% colnames(meta)) {
    message("\n  [E] RNA-seq vs TCR-seq TRBV comparison...")

    # Whole sample
    rna_cmp_whole <- tryCatch(
      compare_rnaseq_vs_tcrseq(meta, "trbv_rna","trbv_gene","sample_id", min_cells),
      error=function(e) NULL
    )

    # Stratified
    rna_cmp_sub <- NULL
    if (has_subset) {
      rna_cmp_sub <- tryCatch(
        compare_rnaseq_vs_tcrseq(
          meta %>% filter(subset %in% ALL_SUBSETS),
          "trbv_rna","trbv_gene", c("sample_id","subset"), 20
        ),
        error=function(e) NULL
      )
    }
    rna_comparison <- list(whole=rna_cmp_whole, by_subset=rna_cmp_sub)

    if (!is.null(rna_cmp_whole)) {
      write.csv(rna_cmp_whole$per_group,
                paste0(out,"_rnaseq_vs_tcrseq_whole.csv"), row.names=FALSE)
      if (nrow(rna_cmp_whole$gene_concordance) > 0)
        write.csv(rna_cmp_whole$gene_concordance,
                  paste0(out,"_gene_concordance_trbv.csv"), row.names=FALSE)

      safe(plot_rnaseq_vs_tcrseq_scatter(rna_cmp_whole,"sample_id",chain_label="TRBV"),
           paste0(out,"_rnaseq_vs_tcrseq_scatter.pdf"), w=7)
      safe(plot_gene_concordance_bar(rna_cmp_whole, chain_label="TRBV"),
           paste0(out,"_gene_concordance_bar.pdf"), w=8, h=7)
    }
    if (!is.null(rna_cmp_sub)) {
      write.csv(rna_cmp_sub$per_group,
                paste0(out,"_rnaseq_vs_tcrseq_by_subset.csv"), row.names=FALSE)
      if (nrow(rna_cmp_sub$per_group) >= 2) {
        p_rna_sub <- ggplot(
          rna_cmp_sub$per_group,
          aes(H_tcrseq, H_rna, color=subset, label=sample_id)
        ) +
          geom_abline(slope=1,intercept=0,color="red",linetype="dotted") +
          geom_point(size=3,alpha=.85) +
          geom_text_repel(size=2.5,max.overlaps=12) +
          facet_wrap(vars(subset),scales="free") +
          scale_color_manual(values=SUBSET_COLORS) +
          labs(title=paste("RNA-seq vs TCR-seq TRBV Usage Diversity by Subset —",dataset_id),
               x="H_TRBV (TCR-seq)", y="H_TRBV (RNA-seq)") +
          theme_bw(base_size=11) + theme(legend.position="none")
        safe(p_rna_sub, paste0(out,"_rnaseq_vs_tcrseq_by_subset.pdf"), w=13, h=6)
      }
    }

    # TRAV comparison
    if ("trav_rna" %in% colnames(meta) && "trav_gene_std" %in% colnames(meta)) {
      rna_cmp_trav <- tryCatch(
        compare_rnaseq_vs_tcrseq(meta, "trav_rna","trav_gene_std","sample_id",min_cells),
        error=function(e) NULL
      )
      if (!is.null(rna_cmp_trav)) {
        write.csv(rna_cmp_trav$per_group,
                  paste0(out,"_rnaseq_vs_tcrseq_trav.csv"), row.names=FALSE)
        safe(plot_gene_concordance_bar(rna_cmp_trav, chain_label="TRAV"),
             paste0(out,"_gene_concordance_trav_bar.pdf"), w=8, h=7)
      }
    }
  }

  # ── [F] TCR expression heatmap ────────────────────────────────────────────
  if (has_expr && !is.null(expr_data$expr)) {
    message("\n  [F] TCR gene expression heatmaps...")

    # Whole-sample heatmap
    tryCatch({
      pdf(paste0(out,"_tcr_expression_heatmap.pdf"), width=12, height=10)
      plot_tcr_expression_heatmap(
        expr_data$expr, meta, group_vars="sample_id",
        chains=c("TRBV","TRAV","TRBJ","TRAJ"),
        title=paste("TCR Gene Expression —",dataset_id,"(all cells)")
      )
      dev.off()
    }, error=function(e) { dev.off(); warning("TCR heatmap: ",e$message) })

    # Per-subset heatmaps
    if (has_subset) {
      for (ss in intersect(ALL_SUBSETS, unique(meta$subset))) {
        sub_m <- meta %>% filter(subset==ss)
        if (nrow(sub_m) < min_cells) next
        tryCatch({
          pdf(paste0(out,"_tcr_expression_heatmap_",ss,".pdf"), width=12, height=10)
          plot_tcr_expression_heatmap(
            expr_data$expr, sub_m, group_vars="sample_id",
            title=paste("TCR Gene Expression —",dataset_id,"—",ss,"T cells")
          )
          dev.off()
        }, error=function(e) { dev.off(); warning("Heatmap ",ss,": ",e$message) })
      }
    }
  }

  # ── Standard plots ────────────────────────────────────────────────────────
  safe(plot_trbv_stacked_bar(meta,"sample_id",20), paste0(out,"_trbv_stacked_bar.pdf"))
  safe(plot_clonotype_rank_abundance(meta,"sample_id"), paste0(out,"_rank_abundance.pdf"),w=9)
  safe(plot_h_trbv_vs_h_clonotype(res_whole$comparison$table,"sample_id"),
       paste0(out,"_h_trbv_vs_h_clonotype.pdf"),w=7)
  safe(plot_clonality_comparison(res_whole$comparison$table,"sample_id"),
       paste0(out,"_clonality_comparison.pdf"),w=7)
  safe(plot_bland_altman(res_whole$comparison$table,"sample_id"),
       paste0(out,"_bland_altman.pdf"),w=7)
  safe(plot_within_trbv_diversity(bm,"sample_id"),
       paste0(out,"_within_trbv_diversity.pdf"),w=14,h=8)

  for (cfg in list(list("trav_gene_std","TRAV"), list("traj_gene_std","TRAJ"),
                   list("trbj_gene","TRBJ"))) {
    if (cfg[[1]] %in% colnames(meta) &&
        sum(!is.na(meta[[cfg[[1]]]]) & meta[[cfg[[1]]]] != "") >= min_cells)
      safe(plot_gene_usage_bar(meta,cfg[[1]],"sample_id",20,
                                title=paste(cfg[[2]],"Usage —",dataset_id)),
           paste0(out,"_",tolower(cfg[[2]]),"_stacked_bar.pdf"))
  }

  if (!is.null(pseudo) && nrow(pseudo$per_level) > 0) {
    sub_var <- if (has_subset) "subset" else NULL
    safe(plot_pseudoclonotype_summary(pseudo, subset_var=sub_var),
         paste0(out,"_pseudoclonotype_summary.pdf"),w=13,h=5)
    safe(plot_pseudo_h_scatter(pseudo,"sample_id",sub_var),
         paste0(out,"_pseudoclonotype_h_scatter.pdf"),w=14,h=6)
  }

  message("  All outputs: ", out, "_*")

  list(
    meta=meta, res_whole=res_whole, sum_whole=sum_whole,
    res_subset=res_subset, sum_subset=sum_subset,
    full_tcr=full_tcr, pseudo=pseudo,
    rna_comparison=rna_comparison, has_expr=has_expr
  )
}


# ── 3. Run all four datasets --------------------------------------------------

blum_tissue_out <- run_dataset_analysis(
  file.path(VAL_DIR,"blum_combined_tissue_data.h5ad"), "blum_tissue")
lechner_out     <- run_dataset_analysis(
  file.path(VAL_DIR,"lechner_thyroiditis_tissue_cd8_seurat_object.h5ad"), "lechner")
luoma_out       <- run_dataset_analysis(
  file.path(VAL_DIR,"luoma_colitis_cd8_seurat_object.h5ad"), "luoma")
blum_pbmc_out   <- run_dataset_analysis(
  file.path(VAL_DIR,"blum_combined_pbmc_data_with_tcr.h5ad"), "blum_pbmc")

all_outputs <- list(
  blum_tissue=blum_tissue_out, lechner=lechner_out,
  luoma=luoma_out,             blum_pbmc=blum_pbmc_out
)


# ── 4. Cross-dataset summaries ------------------------------------------------

message("\n", strrep("=",70), "\nCROSS-DATASET SUMMARY\n", strrep("=",70))

# Whole-sample correlations
corr_whole <- purrr::map_dfr(all_outputs, function(ds) {
  comp <- ds$res_whole$comparison
  if (is.null(comp$spearman)||all(is.na(comp$spearman)))
    return(tibble(rho_H=NA,p_H=NA,rho_clon=NA,n=nrow(comp$table)))
  tibble(rho_H=round(comp$spearman$estimate,3),
         p_H=round(comp$spearman$p.value,4),
         rho_clon=round(comp$clonality_spearman$estimate,3),
         n=nrow(comp$table))
},.id="dataset")
write.csv(corr_whole,
          file.path(OUTPUT_DIR,"cross_dataset_correlation_summary.csv"),row.names=FALSE)
print(corr_whole)

# RNA-seq concordance summary
rna_conc <- purrr::map_dfr(all_outputs, function(ds) {
  if (!ds$has_expr || is.null(ds$rna_comparison$whole)) return(NULL)
  pg <- ds$rna_comparison$whole$per_group
  tibble(mean_concordance=round(mean(pg$concordance,na.rm=TRUE),3),
         n_groups=nrow(pg),
         mean_coverage_rna=round(mean(pg$coverage_rna,na.rm=TRUE),3))
},.id="dataset")
if (nrow(rna_conc)>0) {
  cat("\n--- RNA-seq vs TCR-seq TRBV concordance ---\n")
  print(rna_conc)
  write.csv(rna_conc,
            file.path(OUTPUT_DIR,"cross_dataset_rnaseq_concordance.csv"),row.names=FALSE)
}

# Pseudoclonotype cross-dataset purity
pseudo_all <- purrr::map_dfr(all_outputs, ~ {
  if (is.null(.x$pseudo)||nrow(.x$pseudo$summary)==0) return(NULL)
  .x$pseudo$summary
},.id="dataset")
if (nrow(pseudo_all)>0) {
  write.csv(pseudo_all,
            file.path(OUTPUT_DIR,"cross_dataset_pseudoclonotype_summary.csv"),row.names=FALSE)
  cat("\n--- Cross-dataset pseudoclonotype purity by level ---\n")
  print(pseudo_all %>%
    group_by(pseudo_level) %>%
    summarise(mean_purity=round(mean(mean_purity,na.rm=TRUE),3),
              mean_collision=round(mean(mean_collision,na.rm=TRUE),3),
              datasets=n()) %>%
    arrange(pseudo_level))
}

# Cross-dataset scatter
all_comp <- purrr::map_dfr(all_outputs,~.x$res_whole$comparison$table,.id="dataset")
if (nrow(all_comp)>=3) {
  write.csv(all_comp,
            file.path(OUTPUT_DIR,"cross_dataset_diversity_summary.csv"),row.names=FALSE)
}

# Cross-dataset CD4/CD8/DP/DN scatter
all_subset <- purrr::map_dfr(all_outputs,~ {
  if (is.null(.x$res_subset)) return(NULL)
  .x$res_subset$comparison$table
},.id="dataset")
if (nrow(all_subset)>=3 && "subset" %in% colnames(all_subset)) {
  write.csv(all_subset,
            file.path(OUTPUT_DIR,"cross_dataset_subset_diversity.csv"),row.names=FALSE)
}

message("\nAll outputs saved to: ", OUTPUT_DIR)
