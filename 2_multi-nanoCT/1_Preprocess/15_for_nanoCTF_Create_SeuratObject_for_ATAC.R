#!/usr/bin/env Rscript

# ==============================================================================
# scATAC-seq preprocessing, clustering, annotation, and chromVAR pipeline
# ------------------------------------------------------------------------------
# Purpose:
#   Build a Seurat/Signac object from Cell Ranger ATAC peak counts, perform QC,
#   integrate replicates, cluster cells, annotate clusters using scMultiome-derived
#   cell type marker peaks, run chromVAR, and save the final annotated object.
#
# Workflow:
#   1. Load ATAC peak count matrix and fragment file.
#   2. Create a ChromatinAssay and Seurat object.
#   3. Compute QC metrics and filter cells.
#   4. Run TF-IDF, SVD/LSI, replicate integration, UMAP, and clustering.
#   5. Score cells using cell type-specific marker peaks from scMultiome data.
#   6. Annotate clusters using an external annotation table.
#   7. Add JASPAR motifs, run chromVAR, and save the final object.
#
# Required external files:
#   - filtered_peak_bc_matrix.h5
#   - fragments.tsv.gz
#   - SeuratQCvalues_ATAC.txt
#   - CellTypeAnnotation_ATAC.txt
#   - scMultiome-derived marker peak .rds files
#   - JASPAR2024.sqlite3
#
# Major outputs:
#   - FeaturePlt_CellTypeMarkerPeaks_ATAC_signal_*.png
#   - <REGION>_SO_ATAC.rds
# ==============================================================================

library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(ensembldb)
library(GenomeInfoDb)
library(GenomicRanges)
library(IRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(TFBSTools)
library(chromVAR)
library(motifmatchr)
library(JASPAR2024)
library(dplyr)
library(Matrix)
library(ggplot2)

# config -----------------------------------------------------------------------
REGION <- c("BLA", "Hippo", "Cortex")[1]
PATH_TO_CELLRANGER_AGGR <- "path/to/your/Cellranger/aggr/res"
RESULT_DIR_FINAL <- "/path/to/output/directory"

QC_PARAMETER_FILE <- "path/to/Github/2_multi-nanoCT/1_Preprocess/metadata/SeuratQCvalues_TF_ATAC.txt"
CELLTYPEANNOTATION_FILE <- "path/to/Github/2_multi-nanoCT/1_Preprocess/metadata/CellTypeAnnotation_TF_ATAC.txt"

MULTIOME_DIFFPEAK_DIR <- "path/to/MO/DIFFPEAK" # from 1_scMultiome/2_CellTypeMarkerPeaks/1_CellTypeMarkerPeaks_for_nanoCT_Annotation.R
DIFFPEAK_FILES <- list(
  VGlut  = paste0(MULTIOME_DIFFPEAK_DIR, "/", REGION, "_CellType_MarkerPeaks_Exc_Neu.rds"),
  GABA = paste0(MULTIOME_DIFFPEAK_DIR, "/", REGION, "_CellType_MarkerPeaks_Inh_Neu.rds"),
  Glia  = paste0(MULTIOME_DIFFPEAK_DIR, "/", REGION, "_CellType_MarkerPeaks_Glia.rds")
)

# helper -----------------------------------------------------------------------
Convert_Grange_func <- function(peak){
  
  peak_coords <- do.call(rbind, strsplit(peak, "-"))
  colnames(peak_coords) <- c("chr", "start", "end")
  
  gr_da_peaks <- GRanges(
    seqnames = peak_coords[, "chr"],
    ranges = IRanges(
      start = as.numeric(peak_coords[, "start"]),
      end = as.numeric(peak_coords[, "end"])
    )
  )
  
  return(gr_da_peaks)
  
}

Add_Counts_over_peaks_overlapping_with_CellTypeSpecificPeak <- function(peak_specific, peak_original){
  
  overlap_peak_id <-findOverlaps(peak_original, peak_specific) %>% queryHits()
  h3_counts <- GetAssayData(SO_integrated, assay = "peaks", slot = "data")
  target_counts <- h3_counts[overlap_peak_id, , drop = FALSE]
  lib_size <- Matrix::colSums(SO_integrated, assay = "peaks", slot = "data")
  
  signal_vector <- Matrix::colSums(target_counts) / lib_size
  
  return(signal_vector)
  
}

# load data --------------------------------------------------------------------
count <- Read10X_h5(paste0(PATH_TO_CELLRANGER_AGGR, "/outs/filtered_peak_bc_matrix.h5"))
fragpath <- paste0(PATH_TO_CELLRANGER_AGGR, "/outs/fragments.tsv.gz")

Feature_names <- read.delim(paste0(PATH_TO_CELLRANGER_AGGR, "/outs/filtered_peak_bc_matrix/barcodes.tsv"),
                            header = FALSE,
                            stringsAsFactors = FALSE)

annotation <- GetGRangesFromEnsDb(EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"

# make seurat object -----------------------------------------------------------
chrom_assay_raw <- CreateChromatinAssay(
  counts = count,
  sep = c(":", "-"),
  fragments = fragpath,
  min.cells = 3,
  annotation = annotation
)

SO_raw  <- CreateSeuratObject(
  counts = chrom_assay_raw ,
  assay = "peaks"
)

chrom_assay <- CreateChromatinAssay(
  counts = count,
  sep = c(":", "-"),
  fragments = fragpath,
  min.cells = 3,
  annotation = annotation
)

SO <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)

peaks.keep <- seqnames(granges(SO)) %in% standardChromosomes(granges(SO))
SO <- SO[as.vector(peaks.keep), ]

# perform QC -------------------------------------------------------------------
# compute nucleosome signal score per cell
SO <- NucleosomeSignal(object = SO)

# compute TSS enrichment score per cell
SO <- TSSEnrichment(object = SO)

# filter unwanted cells
qc_params <- read.delim(QC_PARAMETER_FILE, stringsAsFactors = FALSE)
qc_params <- qc_params[which(qc_params$region == REGION),]

SO <- subset(
  x = SO,
  subset = nCount_ATAC < qc_params$nCount_ATAC_max &
    nCount_ATAC > qc_params$nCount_ATAC_min &
    TSS.enrichment > qc_params$TSS.enrichment &
    nucleosome_signal < qc_params$nucleosome_signal)

# Normalization and linear dimensional reduction -------------------------------
SO@meta.data$Rep <- case_when(sub(".+-", "", rownames(SO@meta.data))==1~"Rep1",
                              sub(".+-", "", rownames(SO@meta.data))==2~"Rep2")

SO <- SplitObject(SO, split.by = "Rep")

SO <- lapply(SO, RunTFIDF)
SO <- lapply(SO, FindTopFeatures, min.cutoff = 'q0')
SO <- lapply(SO, RunSVD)

features <- SelectIntegrationFeatures(object.list = SO, nfeatures = 20000)

SO_combined <- merge(SO$Rep2, SO$Rep1)
SO_combined <- FindTopFeatures(SO_combined, min.cutoff = 10)
SO_combined <- RunTFIDF(SO_combined)
SO_combined <- RunSVD(SO_combined)
SO_combined <- RunUMAP(SO_combined, reduction = "lsi", dims = 2:15)

anchors  <- FindIntegrationAnchors(
  object.list = SO,
  anchor.features = features,
  reduction = "rlsi",
  dims = 2:15
)

SO_integrated <- IntegrateEmbeddings(anchorset = anchors, 
                                     new.reduction.name = "integrated.lsi",
                                     reductions = SO_combined[["lsi"]],
                                     dims.to.integrate = 1:15)

SO_integrated <- RunUMAP(SO_integrated, reduction = "integrated.lsi", dims = 2:15)
SO_integrated <- FindNeighbors(object = SO_integrated, reduction = 'integrated.lsi', dims = 2:15)
SO_integrated <- FindClusters(object = SO_integrated, verbose = FALSE, algorithm = 3, resolution = 0.65)

# Annotate clusters using cell type-specific peak ------------------------------
da_peaks_grange_list <- lapply(DIFFPEAK_FILES, 
                               function(dflist){
                                 
                                 lapply(dflist, 
                                        function(df){
                                          tmp <- rownames(df)[which(df$avg_log2F > 1)]
                                          Convert_Grange_func(tmp)
                                        })
                               })

# take overlap with cell type specific peaks
atac_gr <- granges(SO_integrated[["peaks"]])

featurePlt_list <- lapply(da_peaks_grange_list, 
                          function(peaklist){
                            
                            PeakCount_list <- lapply(peaklist, 
                                                     Add_Counts_over_peaks_overlapping_with_CellTypeSpecificPeak, 
                                                     atac_gr)
                            
                            for (i in seq_along(PeakCount_list)){
                              
                              SO_integrated@meta.data[names(PeakCount_list)[i]] <- PeakCount_list[[i]]
                              
                            }
                            
                            FeaturePlot(SO_integrated, features = names(PeakCount_list))
                            
                          })

for(i in 1:3){
  ggsave(paste0(RESULT_DIR_FINAL, "/FeaturePlt_CellTypeMarkerPeaks_ATAC_signal_", c("VGlut", "GABA", "Inh")[i], ".png"),
         featurePlt_list[[i]])
}

# Annotate each cluster --------------------------------------------------------
annotation_params <- read.delim(CELLTYPEANNOTATION_FILE, stringsAsFactors = FALSE)
annotation_params <- annotation_params[which(annotation_params$region == REGION),]

cluster_to_celltype <- setNames(
  annotation_params$CellType,
  annotation_params$ClusterNumber
)

clusters <- as.character(SO_integrated@meta.data$seurat_clusters)

SO_integrated@meta.data$Annotation_first <- cluster_to_celltype[clusters]
SO_integrated@meta.data$Annotation_first <- factor(SO_integrated@meta.data$Annotation_first, levels = annotation_params$CellType)

# Run ChromVAR -----------------------------------------------------------------
pfm <- getMatrixSet(
  x = "JASPAR2024.sqlite3", # https://jaspar.elixir.no/downloads/
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
SO_integrated <- AddMotifs(
  object = SO_integrated,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

SO_integrated <- RunChromVAR(
  object = SO_integrated,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  new.assay.name = "chromvar"
)

saveRDS(SO_integrated, paste0(RESULT_DIR_FINAL, "/", REGION, "_SO_ATAC.rds"))
