#!/usr/bin/env Rscript

# ==============================================================================
# Identify cell-type-enriched ChromVAR motifs in HC nanoCT Seurat objects
# ------------------------------------------------------------------------------
# Purpose:
#   This script identifies ChromVAR motifs enriched in specific neuronal cell
#   types from HC samples and visualizes the top enriched motifs using violin
#   plots across annotated cell types.
#
# Inputs:
#   - Region-specific Seurat object containing all cells:
#       SO_allCells_<REGION>.rds
#   - Cell-type annotations stored in the Seurat object metadata column:
#       Annotation
#   - Simplified cell-type annotations used for differential testing:
#       Annotation_simple
#   - Cell IDs stored in the Seurat object metadata column:
#       CellID
#   - Sample labels stored in the Seurat object metadata column:
#       Sample
#   - ChromVAR scores stored in the Seurat assay:
#       chromvar
#
# Notes:
#   - Edit the CONFIG section before running.
#   - Region-specific neuronal cell types are defined in CellTypes.
#   - The Seurat object is first filtered to HC cells and then to the selected
#     region-specific cell types.
#   - Differential ChromVAR motif enrichment is calculated for each cell type
#     using FindMarkers.
#   - Violin plots are generated for the top enriched motifs in each cell type.
#
# Outputs:
#   - Violin plots of cell-type-enriched ChromVAR scores:
#       ViolinPlt_ChromVARScore_Enriched_in_<CELLTYPE>.png
# ==============================================================================

library(dplyr)
library(purrr)
library(ggplot2)
library(Seurat)
library(Signac)

# config -----------------------------------------------------------------------
PATH_TO_SEURAT_OBJECT <- "path/to/your/SeuratObject"
FINAL_RES_DIR <- "path/to/final/res"
REGION <- c("BLA", "Hippo", "mPFC")[1]

# import data ------------------------------------------------------------------
SO <- readRDS(paste0(PATH_TO_SEURAT_OBJECT, "/SO_allCells_", REGION, ".rds"))
SO_HC <- subset(SO, Sample %in% c("HC-1", "HC-2"))

CellTypes <- if(REGION == "BLA"){
  c("LA_Chst9", "BA", "BLA_Sst")
}else if(REGION == "Hippo"){
  c("DG", "CA1", "CA3")
}else if(REGION == "mPFC"){
  c("L2_3_IT", "L4_5_IT", "L5_ET", "L5_NP", "L6_CT", "L6_IT", "Sst", "Pvalb")
}

# filter Seurat Object ---------------------------------------------------------
SO_selected <- subset(SO_HC, Annotation %in% CellTypes)

# Calculate differential ChromVAR score
Diff_list <- lapply(CellTypes, 
                    function(ct){
                      
                      FindMarkers(SO_selected,
                                  ident.1 = SO_selected@meta.data$CellID[SO_selected@meta.data$Annotation_simple == ct],
                                  assay = "chromvar",
                                  only.pos = T)
                      
                    })

ViolinPlt_list <- lapply(Diff_list, 
                         function(diff_res){
                           
                           VlnPlot(SO_selected, 
                                   features = rownames(diff_res)[1:50], 
                                   group.by = "Annotation")
                           
                         })

for (i in seq_along(ViolinPlt_list)){
  
  ggsave(paste0(FINAL_RES_DIR, "/ViolinPlt_ChromVARScore_Enriched_in_", CellTypes[i], ".png"), ViolinPlt_list[[i]], width = 20, height = 40)
  
}