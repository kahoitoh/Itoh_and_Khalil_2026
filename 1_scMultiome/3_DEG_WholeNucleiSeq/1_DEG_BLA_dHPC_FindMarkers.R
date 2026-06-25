#!/usr/bin/env Rscript

# ==============================================================================
# AP-1 motif activity-based DEG analysis for BLA / hippocampus
# ------------------------------------------------------------------------------
# This script compares RNA expression between AP-1-high and AP-1-low cells using
# chromVAR AP-1 motif scores and Seurat FindMarkers.
#
# Main steps:
#   1. Load annotated Seurat objects.
#   2. Subset target inhibitory or excitatory neuron populations.
#   3. Select AP-1-high and AP-1-low cells within each condition.
#   4. Run FindMarkers between selected cell groups.
#   5. Save DEG tables.
#
# Comparisons:
#   - CondHigh_vs_CondLow
#   - CondHigh_vs_HCLow
#   - HCHigh_vs_HCLow
#
# Notes:
#   - AP-1 motif ID: MA1141.2 by default.
#   - DEG testing uses the SCT assay with the Wilcoxon test.
#   - Analyses are run separately for BLA or hippocampal neuron populations.
# ==============================================================================

library(Seurat)
library(SeuratObject)
library(dplyr)
library(purrr)

# config -----------------------------------------------------------------------
REGION <- c("BLA", "Hippo")[2]
RESULT_DIR_FINAL <- "path/to/final/res"
PATH_to_SO <- "path/to/SO_list.rds"

# helper -----------------------------------------------------------------------
run_FindMarkers_by_AP1_score <- function(
    seurat_obj,
    motif_id = "MA1141.2",
    celltype_column,
    celltypes_to_keep,
    sample_column = "Sample",
    cell_id_column = "CellID",
    percentile_high,
    percentile_low,
    pval_thresh,
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
    mutate(q90 = quantile(ChromVAR_count, percentile_high, na.rm = TRUE)) %>%
    dplyr::filter(ChromVAR_count >= q90)
  
  df_q10 <- df_chromvar_count %>%
    group_by(HC_Cond4h) %>%
    mutate(q10 = quantile(ChromVAR_count, percentile_low, na.rm = TRUE)) %>%
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
    
    FM_res <- FindMarkers(so_trapped,
                          ident.1 = g1,
                          ident.2 = g2,
                          test.use = "wilcox",
                          min.pct = 0.01,
                          assay = "SCT")
    
    FM_res$Class <- case_when(FM_res$avg_log2FC >= 1 & FM_res$p_val_adj < pval_thresh~"Upregulated",
                              FM_res$avg_log2FC <= -1 & FM_res$p_val_adj < pval_thresh~"Downregulated",
                              TRUE~"Unchanged")
    
    FM_res$Gene_name <- rownames(FM_res)
    
    FM_res

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

# Perform FindMarkers ----------------------------------------------------------

if(REGION == "BLA"){
  
  FM_Inh <- run_FindMarkers_by_AP1_score(
    seurat_obj = SO_Inh,
    motif_id = "MA1141.2",
    celltype_column = "Annotation",
    celltypes_to_keep = c("BLA_Sst", "BLA_Vip"), 
    percentile_high = 0.8,
    percentile_low = 0.2,
    pval_thresh = 0.01,
    output_dir = RESULT_DIR_FINAL,
    output_prefix = "BLA_Inh"
  )
  
  FM_Exc <- run_FindMarkers_by_AP1_score(
    seurat_obj = SO_Exc,
    motif_id = "MA1141.2",
    celltype_column = "Annotation",
    celltypes_to_keep = SO_Exc@meta.data$Annotation %>% unique(), 
    percentile_high = 0.9,
    percentile_low = 0.1,
    pval_thresh = 0.001,
    output_dir = RESULT_DIR_FINAL,
    output_prefix = "BLA_Exc"
  )
  
}else if(REGION == "Hippo"){
  
  FM_Inh <- run_FindMarkers_by_AP1_score(
    seurat_obj = SO_Inh,
    motif_id = "MA1141.2",
    celltype_column = "Annotation",
    celltypes_to_keep = c("Sst", "Vip"), 
    percentile_high = 0.8,
    percentile_low = 0.2,
    pval_thresh = 0.01,
    output_dir = RESULT_DIR_FINAL,
    output_prefix = "Hippo_Inh"
  )
  
  FM_DG <- run_FindMarkers_by_AP1_score(
    seurat_obj = SO_Exc,
    motif_id = "MA1141.2",
    celltype_column = "Annotation",
    celltypes_to_keep = "DG", 
    percentile_high = 0.9,
    percentile_low = 0.1,
    pval_thresh = 0.001,
    output_dir = RESULT_DIR_FINAL,
    output_prefix = "Hippo_DG"
  )
  
  FM_CA1 <- run_FindMarkers_by_AP1_score(
    seurat_obj = SO_Exc,
    motif_id = "MA1141.2",
    celltype_column = "Annotation",
    celltypes_to_keep = "CA1", 
    percentile_high = 0.9,
    percentile_low = 0.1,
    pval_thresh = 0.001,
    output_dir = RESULT_DIR_FINAL,
    output_prefix = "Hippo_CA1"
  )
  
  FM_CA3 <- run_FindMarkers_by_AP1_score(
    seurat_obj = SO_Exc,
    motif_id = "MA1141.2",
    celltype_column = "Annotation",
    celltypes_to_keep = "CA3", 
    percentile_high = 0.9,
    percentile_low = 0.1,
    pval_thresh = 0.001,
    output_dir = RESULT_DIR_FINAL,
    output_prefix = "Hippo_CA3"
  )
}