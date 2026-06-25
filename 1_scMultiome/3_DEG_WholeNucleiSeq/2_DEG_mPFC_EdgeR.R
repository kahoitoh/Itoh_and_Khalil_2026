#!/usr/bin/env Rscript

# ==============================================================================
# AP-1 motif activity-based pseudobulk DEG analysis for mPFC
# ------------------------------------------------------------------------------
# This script compares RNA expression between AP-1-high and AP-1-low cells using
# chromVAR motif deviation scores and edgeR pseudobulk analysis.
#
# Main steps:
#   1. Load annotated Seurat objects.
#   2. Subset target inhibitory or excitatory neuron types.
#   3. Select top/bottom 10% AP-1 chromVAR-score cells within each condition.
#   4. Aggregate RNA counts by AP-1 group and replicate.
#   5. Run edgeR and save DEG tables.
#
# Comparisons:
#   - CondHigh_vs_CondLow
#   - CondHigh_vs_HCLow
#   - HCHigh_vs_HCLow
#
# Notes:
#   - AP-1 motif ID: MA1141.2 by default.
#   - Input object list must contain Inh_Neu and Exc_Neu.
#   - Use separate output prefixes for inhibitory and excitatory analyses.
# ==============================================================================

library(dplyr)
library(purrr)
library(Seurat)
library(edgeR)
library(tibble)

# config -----------------------------------------------------------------------
RESULT_DIR_FINAL <- "path/to/Github/1_scMultiome/3_DEG_WholeNucleiSeq/data"
PATH_to_SO <- "/path/to/SeuratObjectList"

# helper -----------------------------------------------------------------------
run_edgeR_by_AP1_score <- function(
    seurat_obj,
    motif_id = "MA1141.2",
    celltype_column,
    celltypes_to_keep,
    sample_column = "Sample",
    cell_id_column = "CellID",
    output_dir,
    output_prefix
) {
  
  seurat_obj$HC_Cond4h <- sub("-[0-9]", "", seurat_obj@meta.data[[sample_column]])
  
  cells_to_keep <- rownames(seurat_obj@meta.data)[
    seurat_obj@meta.data[[celltype_column]] %in% celltypes_to_keep
  ]
  
  seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
  
  chromvar_count_selected <- seurat_obj@assays$chromvar$data[
    rownames(seurat_obj@assays$chromvar$data) == motif_id,
  ]
  
  df_chromvar_count <- data.frame(
    ChromVAR_count = as.numeric(chromvar_count_selected),
    CellID = names(chromvar_count_selected),
    HC_Cond4h = seurat_obj@meta.data$HC_Cond4h
  )
  
  df_q90 <- df_chromvar_count %>%
    group_by(HC_Cond4h) %>%
    mutate(q90 = quantile(ChromVAR_count, 0.9, na.rm = TRUE)) %>%
    dplyr::filter(ChromVAR_count >= q90)
  
  df_q10 <- df_chromvar_count %>%
    group_by(HC_Cond4h) %>%
    mutate(q10 = quantile(ChromVAR_count, 0.1, na.rm = TRUE)) %>%
    dplyr::filter(ChromVAR_count <= q10)
  
  group1_list <- list(
    CondHigh_vs_CondLow = df_q90$CellID[df_q90$HC_Cond4h == "Cond4h"],
    CondHigh_vs_HCLow   = df_q90$CellID[df_q90$HC_Cond4h == "Cond4h"],
    HCHigh_vs_HCLow     = df_q90$CellID[df_q90$HC_Cond4h == "HC"]
  )
  
  group2_list <- list(
    CondHigh_vs_CondLow = df_q10$CellID[df_q10$HC_Cond4h == "Cond4h"],
    CondHigh_vs_HCLow   = df_q10$CellID[df_q10$HC_Cond4h == "HC"],
    HCHigh_vs_HCLow     = df_q10$CellID[df_q10$HC_Cond4h == "HC"]
  )
  
  result_list <- map2(group1_list, group2_list, function(g1, g2) {
    
    so_trapped <- subset(seurat_obj, CellID %in% c(g1, g2))
    
    so_trapped$comparison_name <- case_when(
      so_trapped@meta.data[[cell_id_column]] %in% g1 ~ "g1",
      so_trapped@meta.data[[cell_id_column]] %in% g2 ~ "g2"
    )
    
    so_trapped$comparison_name_with_rep <- paste0(
      so_trapped@meta.data$comparison_name,
      "-",
      sub(".+-", "", so_trapped@meta.data[[sample_column]])
    )
    
    get_group_counts <- function(group_rep) {
      mat <- subset(so_trapped, comparison_name_with_rep == group_rep) %>%
        GetAssayData(assay = "RNA", layer = "counts")
      apply(mat, 1, sum)
    }
    
    counts_pseudobulk <- data.frame(
      group1_1 = get_group_counts("g1-1"),
      group1_2 = get_group_counts("g1-2"),
      group2_1 = get_group_counts("g2-1"),
      group2_2 = get_group_counts("g2-2")
    )
    
    group <- factor(
      c("Group1", "Group1", "Group2", "Group2"),
      levels = c("Group2", "Group1")
    )
    
    edgeR_obj <- DGEList(counts = counts_pseudobulk, group = group)
    keep <- filterByExpr(edgeR_obj)
    edgeR_obj <- edgeR_obj[keep, , keep.lib.sizes = TRUE]
    edgeR_obj <- calcNormFactors(edgeR_obj)
    
    design <- model.matrix(~ group)
    edgeR_obj <- estimateDisp(edgeR_obj, design)
    fit <- glmQLFit(edgeR_obj, design)
    qlf <- glmQLFTest(fit, coef = 2)
    
    topTags(qlf, n = Inf)$table %>%
      as.data.frame() %>%
      rownames_to_column("Gene_name") %>%
      mutate(
        DEG = case_when(
          FDR < 0.05 & logFC > 0 ~ "Upregulated",
          FDR < 0.05 & logFC < 0 ~ "Downregulated",
          TRUE ~ "Unchange"
        )
      )
  })
  
  for (name in names(result_list)) {
    write.table(
      result_list[[name]],
      paste0(output_dir, "/", output_prefix, "_", name, ".txt"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  
  return(result_list)
}

# load data --------------------------------------------------------------------
SO_list <- readRDS(PATH_to_SO)

SO_Inh <- SO_list$Inh_Neu
SO_Exc <- SO_list$Exc_Neu

# Perform EdgeR ----------------------------------------------------------------
EdgeR_Inh <- run_edgeR_by_AP1_score(
  seurat_obj = SO_Inh,
  motif_id = "MA1141.2",
  celltype_column = "Annotation",
  celltypes_to_keep = c("Pvalb", "Sst"),
  output_dir = RESULT_DIR_FINAL,
  output_prefix = "mPFC_Inh"
)

EdgeR_Exc <- run_edgeR_by_AP1_score(
  seurat_obj = SO_Exc,
  motif_id = "MA1141.2",
  celltype_column = "Annotation",
  celltypes_to_keep = SO_Exc@meta.data$Annotation %>% unique(), #c("L2_3_IT", "L4_5_IT", "L5_ET", "L5_NP", "L6_CT", "L6_IT"),
  output_dir = RESULT_DIR_FINAL,
  output_prefix = "mPFC_Exc"
)
