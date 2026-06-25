#!/usr/bin/env Rscript

# ==============================================================================
# Export neuron barcodes from nanoCT Seurat objects
# ------------------------------------------------------------------------------
# Purpose:
#   This script extracts barcodes of annotated neuronal cells from the K27ac
#   nanoCT Seurat object for each brain region. HC1 and HC2 barcodes are exported
#   separately for downstream analyses.
#
# Inputs:
#   - Region-specific combined nanoCT Seurat object:
#       <REGION>_combined_SO_with_celltype_peak_scores.rds
#   - Cell-type annotations stored in the Seurat object metadata column:
#       Annotation
#
# Notes:
#   - Edit the CONFIG section before running.
#   - Neuronal cell types are defined separately for BLA, Hippo, and mPFC.
#   - Barcodes are extracted from the K27ac_ac modality.
#
# Outputs:
#   - <REGION>_Barcode_Neuron_HC1.txt
#   - <REGION>_Barcode_Neuron_HC2.txt
# ==============================================================================

library(Seurat)

# Config -----------------------------------------------------------------------
REGION <- c("BLA", "Hippo", "mPFC")[3]

RESULT_DIR_FINAL <- paste0("path/to/your/res/nCT_", REGION, "_baseres_Cells_Having_All_modalities")

combined_objects <- readRDS(paste0(RESULT_DIR_FINAL, "/", REGION, "_combined_SO_ChromVAR.rds")) # from 2_multi-nanoCT/1_Preprocess/2_for_histone_built_SeuratObject.R

celltype_neuron <- if(REGION == "BLA"){
  c("BA", "LA", "Sst", "Vip")
}else if(REGION == "Hippo"){
  c("CA1", "CA3", "DG", "Inh")
}else if(REGION == "mPFC"){
  c("L2_3_IT", "L4_5_IT", "L5_ET", "L5_NP", "L6_CT", "L6_IT", "Sst_PV", "Vip")
}

# helper -----------------------------------------------------------------------
get_barcode_neuron <- function(SO){
  
  SO <- subset(SO, Annotation %in% celltype_neuron)
  df <- data.frame(barcode = SO@meta.data$barcode,
             sample = sub("[0-9].+", "", rownames(SO@meta.data)))
  df_HC_1 <- df[which(SO@meta.data$orig.ident == "HC1"),]
  df_HC_2 <- df[which(SO@meta.data$orig.ident == "HC2"),]
  return(list(df_HC_1, df_HC_2))
  
}

# get barcode from each replicate ---------------------------------------------
Barcode_Neuron_Cond_HC <- get_barcode_neuron(combined_objects$K27ac_ac)

write.table(Barcode_Neuron_Cond_HC[[1]], paste0(RESULT_DIR_FINAL, "/", REGION, "_Barcode_Neuron_HC1.txt"),
            sep = "\t", quote = F, col.names = F, row.names = F)
write.table(Barcode_Neuron_Cond_HC[[2]], paste0(RESULT_DIR_FINAL, "/", REGION, "_Barcode_Neuron_HC2.txt"),
            sep = "\t", quote = F, col.names = F, row.names = F)
