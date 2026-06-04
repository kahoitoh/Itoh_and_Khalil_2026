#!/usr/bin/env Rscript

# ==============================================================================
# Prepare cell inputs for downstream SCENIC+ analysis
# ------------------------------------------------------------------------------
# This script prepares annotation and RNA count files from annotated multiome
# Seurat objects for downstream regulatory analysis.
# ==============================================================================

library(Seurat)
library(Signac)
library(Matrix)
library(dplyr)

# config -----------------------------------------------------------------------
PATH_TO_OBJECT_LIST <- "PATH/to/SO_list"
FINAL_RES_PATH <- "PATH/to/final/reslt"
REGION <- c("BLA", "Hippo", "mPFC")[1]

CELLTYPE_to_CHOOSE <- if(REGION == "BLA"){
  
  list("BA", "LA_Chst9", "BLA_Sst")
  
}else if(REGION == "Hippo"){
  
  list("DG", "CA1", "CA3")
  
}else if(REGION == "mPFC"){
  
  list(c("L2_3_IT", "L4_5_IT", "L5_ET", "L5_NP", "L6_CT", "L6_IT"),
       "Sst", "Pvalb")
  
}

# helper -----------------------------------------------------------------------
get_anno_table_func <- function(sl){
  
  anno_df_list <- lapply(sl, 
                         function(SO){
                           
                           SO@meta.data$Condition <- sub("-1|-2", "", SO@meta.data$Sample)
                           
                           df <- data.frame(barcode = SO@meta.data$CellID,
                                            Cell = SO@meta.data$Annotation,
                                            Sample = "10x_multiome_brain",
                                            Condition = SO@meta.data$Condition,
                                            CellCondition = paste0(SO@meta.data$Annotation, SO@meta.data$Condition))
                           df
                           
                         })
  
  
  anno_df_All <- rbind(anno_df_list$Exc_Neu,
                       anno_df_list$Inh_Neu,
                       anno_df_list$Glia)
  
  rownames(anno_df_All) <- paste0(anno_df_All$barcode, "-10x_multiome_brain")
  
  anno_df_All
}

# load data --------------------------------------------------------------------
SO_list <- readRDS(paste0(PATH_TO_OBJECT_LIST, "/", REGION, "_SO_list_Annotated.rds"))

# Prepare annotation file containing whole cells for pycisTopic consensus peak calling
anno_df_whole <- get_anno_table_func(SO_list)
write.table(anno_df_whole, paste0(FINAL_RES_PATH, "/Whole_", REGION, "_anno.tsv"), sep = "\t", quote = F)

# Prepare TRAPed Seurat object
SO_TRAP_list <- lapply(SO_list, 
                       function(SO){
                         
                         SO_Cond4h <- subset(SO, subset = Sample %in% c("Cond4h-1", "Cond4h-2"))
                         
                         chromvar_count_selected <- SO_Cond4h@assays$chromvar$data[which(rownames(SO_Cond4h@assays$chromvar$data) == "MA1141.2"),]
                         df_chromvar_count <- data.frame(ChromVAR_count = chromvar_count_selected)
                         df_chromvar_count$CellID <- rownames(df_chromvar_count)
                         df_chromvar_count$Celltype <- SO_Cond4h@meta.data$Annotation
                         
                         df_chromvar_count_High <- df_chromvar_count %>% 
                           group_by(., Celltype) %>% 
                           mutate(q80 = quantile(ChromVAR_count, 0.8, na.rm = TRUE)) %>%
                           dplyr::filter(ChromVAR_count >= q80)
                         
                         df_chromvar_count_High$Condition <- "High"
                         
                         df_chromvar_count_Low <- df_chromvar_count %>% 
                           group_by(., Celltype) %>% 
                           mutate(q20 = quantile(ChromVAR_count, 0.2, na.rm = TRUE)) %>%
                           dplyr::filter(ChromVAR_count <= q20)
                         
                         df_chromvar_count_Low$Condition <- "Low"
                         
                         SO_Cond_High_Low <- subset(SO, CellID %in% c(df_chromvar_count_Low$CellID, df_chromvar_count_High$CellID))
                         SO_Cond_High_Low@meta.data$High_low <- case_when(SO_Cond_High_Low@meta.data$CellID %in% df_chromvar_count_Low$CellID~"Low",
                                                                          TRUE~"High")
                         
                         SO_Cond_High_Low
                         
                       })

# Prepare annotation file containing TRAPed each cell type for constracting PycisTopic Object
anno_df_TRAP <- get_anno_table_func(SO_TRAP_list)

anno_df_TRAP_celltype_list <- lapply(CELLTYPE_to_CHOOSE, 
                                     function(ct){
                                       
                                       anno_df_TRAP[which(anno_df_TRAP$Cell %in% ct),]
                                       
                                     }, anno_df_TRAP)

# Prepare RNA raw counts containing TRAPed each cell type for SCENIC+
RNACount_METAdata_TRAP_celltype_list <- lapply(CELLTYPE_to_CHOOSE, 
                                      function(ct){
                                        
                                        if(ct %in% c("BLA_Sst", "Sst", "Pvalb")){
                                          
                                          tmp_SO <- SO_TRAP_list$Inh_Neu
                                          
                                        }else{
                                          tmp_SO <- SO_TRAP_list$Exc_Neu
                                        }
                                        
                                        # get RNA count
                                        CellID_ct <- tmp_SO$CellID[which(tmp_SO$Annotation %in% ct)]
                                        
                                        mat_all <- GetAssayData(tmp_SO, assay = "RNA", slot = "counts")
                                        
                                        mat_selected <- mat_all[,CellID_ct]
                                        
                                        # get metadata
                                        meta <- tmp_SO@meta.data
                                        meta$barcode <- rownames(meta)
                                        meta$Cells_condion <- paste0(meta$Annotation, "-", meta$High_low)
                                        meta$Cells <- meta$Annotation
                                        
                                        res_list <- list(mat_selected, meta)
                                        res_list
                                        
                                      }, SO_TRAP_list)

# Export -----------------------------------------------------------------------
for (i in seq_along(CELLTYPE_to_CHOOSE)){
  write.table(anno_df_TRAP_celltype_list[[i]], paste0(FINAL_RES_PATH, "/TRAPed_", REGION, "_anno.tsv"), sep = "\t", quote = F)
  
  writeMM(RNACount_METAdata_TRAP_celltype_list[[i]][[1]], file = paste0(FINAL_RES_PATH, "/TRAPed_", REGION, "_", CELLTYPE_to_CHOOSE[i], "_rna_counts.mtx"))
  
  write.table(colnames(RNACount_METAdata_TRAP_celltype_list[[i]][[1]]), file = paste0(FINAL_RES_PATH, "TRAPed_", REGION, "_", CELLTYPE_to_CHOOSE[i], "_barcodes.tsv"), 
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  write.table(rownames(RNACount_METAdata_TRAP_celltype_list[[i]][[1]]), file = paste0(FINAL_RES_PATH, "TRAPed_", REGION, "_", CELLTYPE_to_CHOOSE[i], "_genes.tsv"), 
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  write.table(rownames(RNACount_METAdata_TRAP_celltype_list[[i]][[2]]), file = paste0(FINAL_RES_PATH, "TRAPed_", REGION, "_", CELLTYPE_to_CHOOSE[i], "_metadata.tsv"), 
              quote = FALSE, sep = "\t", row.names = FALSE)
}
