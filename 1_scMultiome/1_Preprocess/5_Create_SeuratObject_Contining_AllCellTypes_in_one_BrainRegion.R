#!/usr/bin/env Rscript

# ==============================================================================
# Generate peak assay and chromVAR motif activity scores for all singlet cells
# ------------------------------------------------------------------------------
# Purpose:
#   This script prepares a combined Seurat object containing RNA, ATAC peak counts,
#   motif annotations, and chromVAR motif deviation scores for downstream analyses.
#
# Main steps:
#   1. Load annotated Seurat object lists and the original QC-filtered Seurat object.
#   2. Retain singlet cells assigned to excitatory neurons, inhibitory neurons,
#      or glial cells.
#   3. Call ATAC peaks using MACS2 on the filtered cell set.
#   4. Remove peaks on non-standard chromosomes and mm10 blacklist regions.
#   5. Quantify peak accessibility and add a peak assay to the Seurat object.
#   6. Normalize RNA counts using SCTransform.
#   7. Add JASPAR motif annotations to the peak assay.
#   8. Calculate chromVAR motif deviation scores.
#   9. Save the processed all-cell Seurat object.
#
# Inputs:
#   - <REGION>_SO_list_Annotated.rds
#   - <REGION>_SO_SeuratQC.rds
#   - 10x Multiome ATAC fragment file
#   - JASPAR2024 motif database
#
# Outputs:
#   - SO_allCells_<REGION>.rds
#
# Notes:
#   - Genome build: mm10.
#   - Peak calling is performed using MACS2.
#   - Motif annotations are based on JASPAR2024 vertebrate CORE motifs.
#   - chromVAR scores are stored in the chromvar assay.
# ==============================================================================

library(dplyr)
library(purrr)
library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)

# config
REGION <- c("BLA", "Hippo", "mPFC")[1]
ORIGINAL_SeuratObjectList_DIR <- "oath/to/your/SO_list"
FINAL_RES_DIR <- "path/to/your/final/result"

# import Seurat object list ----------------------------------------------------
SO_list <- readRDS(paste0(ORIGINAL_SeuratObjectList_DIR, "/", REGION, "_SO_list_Annotated.rds"))

CellID_list <- lapply(SO_list,
                      function(SO){
                        
                        SO@meta.data$CellID
                        
                      })

SO_SeuratQC <- readRDS(paste0(ORIGINAL_SeuratObjectList_DIR, "/", REGION, "_SO_SeuratQC.rds"))
SO_SeuratQC@meta.data$CellID <- rownames(SO_SeuratQC@meta.data)

# get singlet ------------------------------------------------------------------
SO_all <- subset(SO_SeuratQC, CellID %in% c(CellID_list$Exc_Neu, CellID_list$Inh_Neu, CellID_list$Glia))

# peak calling -----------------------------------------------------------------
DefaultAssay(SO_all) <- "ATAC"
peaks <- CallPeaks(SO_all, macs2.path = "/path/to/your/macs2")

annotation <- GetGRangesFromEnsDb(EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"

fragpath <- "path/to/your/atac_fragments.tsv.gz"

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = Signac::blacklist_mm10, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(SO_all),
  features = peaks,
  cells = colnames(SO_all)
)

SO_all[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

DefaultAssay(SO_all) <- "peaks"

# Normalize RNA ----------------------------------------------------------------
if(REGION == "mPFC"){
  SO_all@assays$RNA <- split(x = SO_all@assays$RNA, f = SO_all@meta.data$Sample)
}

DefaultAssay(SO_all) <- "RNA"

SO_all <- SCTransform(SO_all)


# Calculate ChromVAR score -----------------------------------------------------
pfm <- getMatrixSet(
  x = "JASPAR2024.sqlite3", # https://jaspar.elixir.no/downloads/
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

SO_all <- AddMotifs(
  object = SO_all,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

SO_all <- RunChromVAR(
  object = SO_all,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  new.assay.name = "chromvar"
)

saveRDS(SO_all, paste0(FINAL_RES_DIR, "/SO_allCells_", REGION, ".rds"))
