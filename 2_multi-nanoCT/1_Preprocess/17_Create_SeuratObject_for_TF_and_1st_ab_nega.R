#!/usr/bin/env Rscript

# ==============================================================================
# TF and 1st antibody negative ATAC object construction
# ------------------------------------------------------------------------------
# Build Seurat/Signac objects from TF and 1st antibody negative fragment files,
# quantify reads over broadPeak regions, keep common cells between the two
# samples, add the 1st antibody negative assay to the TF object, and transfer
# cell-type annotation and FOS motif scores from the annotated scATAC object.
#
# Main steps:
#   1. Load broadPeak regions and annotated scATAC object
#   2. Create peak count matrices from TF and 1st antibody negative fragments
#   3. Keep cells shared between TF and 1st antibody negative datasets
#   4. Add 1st antibody negative peak assay to the TF object
#   5. Transfer annotation, replicate ID, and FOS chromVAR motif score
#   6. Save annotated TF / 1st antibody negative Seurat object
#
# Inputs:
#   - TF fragment file
#   - 1st antibody negative fragment file
#   - broadPeak file
#   - annotated scATAC Seurat object with chromVAR scores
#
# Output:
#   - SO_REGION_TF_1st_ab_nega_Annotated.rds
# ==============================================================================

library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)
library(dplyr)
library(Matrix)

# config -----------------------------------------------------------------------
REGION <- c("BLA", "Hippo", "Cortex")[1]

fragments_TF <- paste0("../../Fragments_", REGION, "_TF.sorted.tsv.gz")   
fragments_1st_ab_nega <- paste0("../../Fragments_", REGION, "_1st_ab_nega.sorted.tsv.gz") 

BROADPEAK_FILE <- "path/to/broadPeak/NA_peaks.broadPeak"  
PATH_TO_SO_ATAC <- "path/to/PATH/TO/SO/ATAC"
RESULT_DIR_FINAL <- "/path/to/output/directory"
# helper -----------------------------------------------------------------------
Make_SO <- function(fragments, peaks_df){
  
  peaks_gr <- makeGRangesFromDataFrame(peaks_df, 
                                       keep.extra.columns = TRUE,
                                       seqnames.field = "chr",
                                       start.field = "start",
                                       end.field   = "end")
  
  tmp <- CreateFragmentObject(
    fragments,
    verbose = TRUE)

  peak_mat <- FeatureMatrix(
    fragments = tmp,
    features  = peaks_gr
  )
  
  chrom_assay <- CreateChromatinAssay(
    counts    = peak_mat,
    fragments = fragments,
    sep       = c(":", "-")
  )
  
  CreateSeuratObject(counts = chrom_assay, assay = "peaks", project = "ATAC")
  
}

# import data ------------------------------------------------------------------
peaks_df <- read.table(BROADPEAK_FILE, header = FALSE)
colnames(peaks_df)[1:3] <- c("chr","start","end")

SO_ATAC <- readRDS(paste0(PATH_TO_SO_ATAC, "/", REGION, "_SO_ATAC.rds"))
SO_ATAC@meta.data$Annotation_clean <- sub("_1|_2", "", SO_ATAC@meta.data$Annotation)

# create seurat object for TF and 1st ab nega ----------------------------------
SO_TF <- Make_SO(fragments_TF, peaks_df)
SO_1st_ab_nega <- Make_SO(fragments_1st_ab_nega, peaks_df)

SO_TF@meta.data$CellID <- rownames(SO_TF@meta.data)
SO_1st_ab_nega@meta.data$CellID <- rownames(SO_1st_ab_nega@meta.data)

Cell_common <- rownames(SO_1st_ab_nega@meta.data)[which(rownames(SO_1st_ab_nega@meta.data) %in% rownames(SO_TF@meta.data))]

SO_TF <- subset(SO_TF, CellID %in% Cell_common)
SO_1st_ab_nega <- subset(SO_1st_ab_nega, CellID %in% Cell_common)

SO_TF[["1st_ab_nega"]] <- SO_1st_ab_nega@assays$peaks

# Add Annotation ---------------------------------------------------------------
df_anno_2 <- data.frame(CellID = SO_ATAC@meta.data$CellID,
                        Annotation = SO_ATAC@meta.data$Annotation_clean,
                        FOS_MotifScore = SO_ATAC@assays$chromvar@data["MA1141.2",])
df_anno_2$CellID <- sub("-2|-1", "", df_anno_2$CellID)
df_anno_2$Rep <- sub(".+-", "", SO_ATAC@meta.data$CellID)

CellID_dups <- df_anno_2$CellID[which(duplicated(df_anno_2$CellID))]

df_anno_2 <- df_anno_2[-which(df_anno_2$CellID %in% CellID_dups),]

SO_TF <- subset(SO_TF, CellID %in% df_anno_2$CellID)

df_anno <- data.frame(CellID = SO_TF@meta.data$CellID)
df_anno <- left_join(df_anno, df_anno_2)

SO_TF <- subset(SO_TF, CellID %in% df_anno$CellID)

SO_TF@meta.data$Annotation <- df_anno$Annotation
SO_TF@meta.data$Rep <- df_anno$Rep
SO_TF@meta.data$FOS_MotifScore <- df_anno$FOS_MotifScore

saveRDS(SO_TF, paste0(RESULT_DIR_FINAL, "/SO_", REGION, "_TF_1st_ab_nega_Annotated.rds"))
