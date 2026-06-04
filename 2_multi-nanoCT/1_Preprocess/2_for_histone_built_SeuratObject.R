#!/usr/bin/env Rscript

# ==============================================================================
# nanoCT ATAC / histone-mark chromatin analysis pipeline
# ------------------------------------------------------------------------------
# Purpose:
#   Build modality-specific Seurat objects from nanoCT ATAC/H3K27ac/H3K27me3
#   fragment files, perform QC, clustering, annotation using overlap analysis 
#   with cell-type-specific peaks, and ChromVAR.
#
# Notes for reproducibility:
#   - Edit only the CONFIG section before running on a new machine.
#   - Sample and modality order are fixed explicitly to avoid order changes from
#     dir()/list.files().
#   - HC1 is treated as a sample without ATAC_aa, matching the original analysis.
# ===============================================================================

# ---- 0. CONFIG ----------------------------------------------------------------

BASE_DIR <- "D:/in_vivo_nanoCT_2502/nCT_Hippo_basedata_rm_HC1_ATAC"
SCRIPT_DIR <- "D:/in_vivo_nanoCT_2502/script"

RESULT_DIR_BROAD <- file.path("..", "nCT_Hippo_baseres_ATACbroadPeaks")
RESULT_DIR_FINAL <- file.path("..", "nCT_Hippo_baseres_rm_HC1_ATAC_Cells_Having_All_modalities")

GENOME <- "mm10"
ASSAY_NAME <- "peaks"
MIN_CELLS_PER_FEATURE <- 1
MIN_FEATURES_PER_CELL <- 1

SAMPLES <- c("Cond1", "Cond2", "HC1", "HC2")
MODALITIES <- c("ATAC_aa", "K27ac_ac", "K27me3_me")
SAMPLES_WITHOUT_ATAC <- c("HC1")

QC_QUANTILE_HIGH <- 0.99
QC_QUANTILE_LOW <- 0.01
FEATURE_MATRIX_PROCESS_N <- 20000

MULTIOME_DIFFPEAK_DIR <- "D:/in_vivo_Multiome_EngramProject/BaseAnalysis_Hippo/5_ATAC_comparison"
DIFFPEAK_FILES <- list(
  GABA  = file.path(MULTIOME_DIFFPEAK_DIR, "GABA_da_peaks_celltype_new_260420.rds"),
  VGlut = file.path(MULTIOME_DIFFPEAK_DIR, "VGlut_da_peaks_celltype_new_260420.rds"),
  Glia  = file.path(MULTIOME_DIFFPEAK_DIR, "non_Neuro_da_peaks_celltype.rds")
)

setwd(BASE_DIR)
dir.create(RESULT_DIR_BROAD, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULT_DIR_FINAL, recursive = TRUE, showWarnings = FALSE)

PARAMETER_FILE <- "D:/in_vivo_Multiome_EngramProject/Github/2_multi-nanoCT/1_Preprocess/config/Parameters_for_clustering_histone.txt"
MARKERPEAKS_THRESHOLD_FILE <- "D:/in_vivo_Multiome_EngramProject/Github/2_multi-nanoCT/1_Preprocess/config/celltype_marker_peak_thresholds.txt"
CELLTYPEANNOTATION_FILE <- "D:/in_vivo_Multiome_EngramProject/Github/2_multi-nanoCT/1_Preprocess/config/CelltypeAnnotation.txt"
REGION <- "Hippo"
# ---- 1. PACKAGES AND CUSTOM FUNCTIONS ----------------------------------------

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  library(future)
  library(stringr)
  library(ggplot2)
  library(gghalves)
  library(ggpubr)
  library(EnsDb.Mmusculus.v79)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(ComplexUpset)
  library(regioneR)
  library(scales)
  library(ggVennDiagram)
  library(dplyr)
  library(Matrix)
})

source(file.path(SCRIPT_DIR, "functions_scCT.R"))

GENOME_ANN <- EnsDb.Mmusculus.v79
BS_GENOME <- BSgenome.Mmusculus.UCSC.mm10

# ---- 2. HELPER FUNCTIONS ------------------------------------------------------

experiment_id <- function(modality, sample) {
  paste0(modality, "_", sample)
}

modalities_for_sample <- function(sample) {
  if (sample %in% SAMPLES_WITHOUT_ATAC) {
    setdiff(MODALITIES, "ATAC_aa")
  } else {
    MODALITIES
  }
}

modality_prefix <- function(modality) {
  strsplit(modality, "_")[[1]][1]
}

read_peak_table <- function(sample, modality) {
  # Original analysis used a dummy ATAC peak for HC1 because HC1 lacks ATAC data.
  # This branch is retained for compatibility, but HC1 ATAC is excluded downstream.
  if (sample == "HC1" && modality == "ATAC_aa") {
    return(data.frame(chr = "chr1", start = 3119747, end = 3120607))
  }
  
  peak_file <- file.path(
    sample,
    modality,
    "peaks",
    "macs_broad",
    paste0(modality_prefix(modality), "_peaks.broadPeak")
  )
  
  peaks <- read.table(peak_file)[, 1:3]
  colnames(peaks) <- c("chr", "start", "end")
  peaks
}

read_metadata <- function(sample, modality, filename = "metadata.csv") {
  metadata_file <- file.path(sample, modality, "cell_picking", filename)
  metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
  rownames(metadata) <- metadata$barcode
  metadata
}


get_modality_from_experiment_id <- function(experiment) {
  # Experiment IDs are formatted like ATAC_aa_Cond1, K27ac_ac_HC2, etc.
  paste(str_split_fixed(experiment, "_", 3)[1, 1:2], collapse = "_")
}

get_sample_from_experiment_id <- function(experiment) {
  str_split_fixed(experiment, "_", 3)[1, 3]
}

subset_to_common_cells <- function(objects, sample) {
  sample_modalities <- modalities_for_sample(sample)
  cell_sets <- lapply(sample_modalities, function(mod) {
    Cells(objects[[experiment_id(mod, sample)]])
  })
  
  common_cells <- Reduce(intersect, cell_sets)
  
  subsetted <- list()
  for (mod in sample_modalities) {
    id <- experiment_id(mod, sample)
    subsetted[[id]] <- subset(objects[[id]], cells = common_cells)
  }
  
  subsetted
}

merge_modality_objects <- function(objects, modality) {
  sample_names <- SAMPLES[!SAMPLES %in% if (modality == "ATAC_aa") SAMPLES_WITHOUT_ATAC else character(0)]
  object_names <- experiment_id(modality, sample_names)
  
  merge(
    x = objects[[object_names[1]]],
    y = objects[object_names[-1]],
    add.cell.ids = sample_names
  )
}

plot_modality_umap <- function(object, title, group.by = NULL, label = FALSE, alpha = 1) {
  p <- DimPlot(
    object,
    group.by = group.by,
    label = label,
    alpha = alpha
  ) +
    ggtitle(title) +
    theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))
  
  if (label) p <- p + NoLegend()
  p
}

add_counts_overlapping_celltype_peaks <- function(combined_objects, modality, original_peaks, specific_peaks) {
  overlap_peak_id <- findOverlaps(original_peaks, specific_peaks) %>% queryHits()
  
  assay_data <- GetAssayData(combined_objects[[modality]], assay = ASSAY_NAME, slot = "data")
  target_counts <- assay_data[overlap_peak_id, , drop = FALSE]
  lib_size <- Matrix::colSums(assay_data)
  
  Matrix::colSums(target_counts) / lib_size
}

parse_dims <- function(x) {
  parts <- strsplit(x, ":")[[1]]
  as.integer(parts[1]):as.integer(parts[2])
}

get_param <- function(region, modality, column) {
  hit <- analysis_params[
    analysis_params$region == region &
      analysis_params$modality == modality,
  ]
  
  if (nrow(hit) != 1) {
    stop("Expected one parameter row for region=", region,
         ", modality=", modality, ", found ", nrow(hit))
  }
  
  hit[[column]]
}

get_marker_threshold <- function(region, celltype, thresholds) {
  hit <- thresholds[
    thresholds$region == region &
      thresholds$celltype == celltype,
    ,
    drop = FALSE
  ]
  
  if (nrow(hit) != 1) {
    stop(
      "Expected one threshold row for region=", region,
      ", celltype=", celltype,
      ", found ", nrow(hit)
    )
  }
  
  hit
}


select_marker_peak_ids <- function(marker_list, region, celltype, thresholds) {
  threshold <- get_marker_threshold(region, celltype, thresholds)
  
  avg_log2fc_min <- as.numeric(threshold$avg_log2FC_min)
  
  lapply(marker_list, function(df) {
    df <- as.data.frame(df)
    
    keep <- df$avg_log2FC > avg_log2fc_min
    
    rownames(df)[keep]
  })
}


marker_id_to_granges <- function(peak_id) {
  
  coords <- do.call(rbind, strsplit(peak_id, "-"))
  
  gr <- GRanges(
    seqnames = coords[, 1],
    ranges = IRanges(
      start = as.integer(coords[, 2]),
      end   = as.integer(coords[, 3])
    )
  )
  
  gr$peak_id <- peak_id
  
  gr
}

# ---- 3. LOAD PEAKS AND CREATE COMMON PEAK SETS --------------------------------

input_peaks <- list()

for (sample in SAMPLES) {
  cat("Loading peaks for sample", sample, "\n")
  
  for (modality in modalities_for_sample(sample)) {
    cat("\t", modality, "(MACS broad peaks)\n")
    input_peaks[[experiment_id(modality, sample)]] <- read_peak_table(sample, modality)
  }
}

input_granges <- lapply(input_peaks, makeGRangesFromDataFrame)

combined_peaks <- list()
for (modality in MODALITIES) {
  modality_samples <- SAMPLES[vapply(
    SAMPLES,
    function(sample) modality %in% modalities_for_sample(sample),
    logical(1)
  )]
  modality_ids <- experiment_id(modality, modality_samples)
  modality_ids <- modality_ids[modality_ids %in% names(input_granges)]
  
  if (length(modality_ids) == 0) {
    combined_peaks[[modality]] <- GenomicRanges::GRanges()
    next
  }
  
  gr <- unlist(
    GenomicRanges::GRangesList(input_granges[modality_ids]),
    use.names = FALSE
  )
  
  combined_peaks[[modality]] <- GenomicRanges::reduce(gr)
  
  seqlevelsStyle(combined_peaks[[modality]]) <- "UCSC"
  
  combined_peaks[[modality]] <- keepStandardChromosomes(
    combined_peaks[[modality]],
    pruning.mode = "coarse"
  )
}

# ---- 4. LOAD METADATA AND CREATE FRAGMENT OBJECTS -----------------------------

metadata_list <- list()

for (sample in SAMPLES) {
  cat("Loading metadata for sample", sample, "\n")
  
  for (modality in modalities_for_sample(sample)) {
    cat("\t", modality, "\n")
    metadata_list[[experiment_id(modality, sample)]] <- read_metadata(sample, modality)
  }
}

metadata_list <- lapply(metadata_list, function(x){x[x$passedMB,]})

fragment_list <- list()

for (sample in SAMPLES) {
  cat("Creating fragment objects for sample", sample, "\n")
  
  for (modality in modalities_for_sample(sample)) {
    id <- experiment_id(modality, sample)
    cat("\t", modality, "\n")
    
    fragment_list[[id]] <- CreateFragmentObject(
      path = file.path(sample, modality, "cellranger", "outs", "fragments.tsv.gz"),
      cells = metadata_list[[id]]$barcode
    )
  }
}

# ---- 5. QUANTIFY PEAKS AND CREATE SEURAT OBJECTS ------------------------------

counts_list <- list()

for (experiment in names(fragment_list)) {
  modality <- get_modality_from_experiment_id(experiment)
  cat("Analysing experiment:", experiment, "\n")
  cat("\tPeaks from modality:", modality, "\n")
  
  counts_list[[experiment]] <- FeatureMatrix(
    fragments = fragment_list[[experiment]],
    features = combined_peaks[[modality]],
    cells = metadata_list[[experiment]]$barcode,
    process_n = FEATURE_MATRIX_PROCESS_N
  )
}

object_list <- list()

for (experiment in names(counts_list)) {
  sample <- get_sample_from_experiment_id(experiment)
  modality <- get_modality_from_experiment_id(experiment)
  
  chrom_assay <- CreateChromatinAssay(
    counts = counts_list[[experiment]],
    fragments = fragment_list[[experiment]],
    genome = GENOME,
    min.cells = MIN_CELLS_PER_FEATURE,
    min.features = MIN_FEATURES_PER_CELL
  )
  
  object_list[[experiment]] <- CreateSeuratObject(
    counts = chrom_assay,
    assay = ASSAY_NAME,
    meta.data = metadata_list[[experiment]],
    project = sample
  )
  
  object_list[[experiment]]$dataset <- experiment
  object_list[[experiment]]$modality <- modality
  object_list[[experiment]]$sample <- sample
}

# ---- 6. QUALITY CONTROL -------------------------------------------------------

object_list_qc <- list()

for (experiment in names(object_list)) {
  logumi_high <- quantile(object_list[[experiment]]$logUMI, QC_QUANTILE_HIGH)
  logumi_low <- quantile(object_list[[experiment]]$logUMI, QC_QUANTILE_LOW)
  peak_ratio_low <- quantile(object_list[[experiment]]$peak_ratio_MB, QC_QUANTILE_LOW)
  
  object_list_qc[[experiment]] <- subset(
    object_list[[experiment]],
    logUMI > logumi_low &
      logUMI < logumi_high &
      peak_ratio_MB > peak_ratio_low
  )
  
  old_n_cell <- nrow(object_list[[experiment]][[]])
  new_n_cell <- nrow(object_list_qc[[experiment]][[]])
  discarded <- old_n_cell - new_n_cell
  
  cat(experiment, "\n")
  cat("\tdiscarded", discarded, "cells (", round(discarded / old_n_cell * 100, 2), "%)\n")
}


# ---- 7. SUBSET TO CELLS PRESENT IN ALL AVAILABLE MODALITIES -------------------

object_list_qc_subset <- list()
for (sample in SAMPLES) {
  object_list_qc_subset <- c(object_list_qc_subset, subset_to_common_cells(object_list_qc, sample))
}

# ---- 8. MERGE SEURAT OBJECTS BY MODALITY --------------------------------------

combined_objects <- list()
for (modality in MODALITIES) {
  combined_objects[[modality]] <- merge_modality_objects(object_list_qc_subset, modality)
}

# ---- 9. NORMALIZATION AND DIMENSIONAL REDUCTION ------------------------------

combined_objects <- lapply(combined_objects, RunTFIDF)
# Warning "Some features contain 0 total counts" can appear after QC filtering.
# This is expected and should be safe to proceed.

combined_objects <- lapply(combined_objects, FindTopFeatures)
combined_objects <- lapply(combined_objects, RunSVD)


# ---- 10. OVERLAP WITH CELL-TYPE-SPECIFIC PEAKS --------------------------------
atac_gr <- granges(combined_objects$ATAC_aa[[ASSAY_NAME]])
k27ac_gr <- granges(combined_objects$K27ac_ac[[ASSAY_NAME]])
k27me3_gr <- granges(combined_objects$K27me3_me[[ASSAY_NAME]])

celltype_peak <- lapply(DIFFPEAK_FILES, readRDS)
MARKER_THRESHOLD_FILE <- read.table(MARKERPEAKS_THRESHOLD_FILE, header = T)


selected_peak_ids <- list()

for (celltype in names(celltype_peak)) {
  
  selected_peak_ids[[celltype]] <- select_marker_peak_ids(
    marker_list = celltype_peak[[celltype]],
    region = REGION,
    celltype = celltype,
    thresholds = MARKER_THRESHOLD_FILE
  )
}



celltype_peak_granges <- lapply(selected_peak_ids,
                                function(list){
                                  
                                  lapply(list, marker_id_to_granges)
                                  
                                })

peak_granges_by_modality <- list(
  ATAC_aa = atac_gr,
  K27ac_ac = k27ac_gr,
  K27me3_me = k27me3_gr
)


for (celltype in names(celltype_peak_granges)) {
  for(subcelltype in names(celltype_peak_granges[[celltype]])){
    
    metadata_column <- paste0(celltype, "_", subcelltype, "Peaks_Counts")
    
    combined_objects$K27ac_ac@meta.data[[metadata_column]] <- add_counts_overlapping_celltype_peaks(
      combined_objects = combined_objects,
      modality = "K27ac_ac",
      original_peaks = peak_granges_by_modality$K27ac_ac,
      specific_peaks = celltype_peak_granges[[celltype]][[subcelltype]]
    )
    
  }
}

# ---- 11. UMAP AND CLUSTERING -----------------------------------------
analysis_params <- read.delim(PARAMETER_FILE, stringsAsFactors = FALSE)

for (modality in MODALITIES) {
  dims <- parse_dims(get_param(REGION, modality, "lsi_dims"))
  resolution <- as.numeric(get_param(REGION, modality, "cluster_resolution"))
  
  combined_objects[[modality]] <- RunUMAP(
    combined_objects[[modality]],
    reduction = "lsi",
    dims = dims
  )
  
  combined_objects[[modality]] <- FindNeighbors(
    combined_objects[[modality]],
    reduction = "lsi",
    dims = dims
  )
  
  combined_objects[[modality]] <- FindClusters(
    combined_objects[[modality]],
    verbose = FALSE,
    algorithm = 3,
    resolution = resolution
  )
}

# ---- 12.Visualize marker peaks H3K27ac signal --------------------------------
lapply(c("VGlut", "GABA", "Glia"), 
       function(ct){
         
         g <- FeaturePlot(combined_objects$K27ac_ac, 
                          features = colnames(combined_objects$K27ac_ac@meta.data)[grep(ct, colnames(combined_objects$K27ac_ac@meta.data))])
         ggsave(paste0(RESULT_DIR_FINAL, "/featurePlt_MarkerPeaks_", ct, ".png"),
                g, width = 10, height = 10)
         
       })

ggsave(paste0(RESULT_DIR_FINAL, "/Dimplt_Clustering.png"),
       DimPlot(combined_objects$K27ac_ac, label = T), width = 4, height = 4)

# ---- 13. Annotate clusters ---------------------------------------------------
annotation_params <- read.delim(CELLTYPEANNOTATION_FILE, stringsAsFactors = FALSE)
annotation_params <- annotation_params[which(annotation_params$region == REGION),]

cluster_to_celltype <- setNames(
  annotation_params$CellType,
  annotation_params$ClusterNumber
)

clusters <- as.character(combined_objects$K27ac_ac@meta.data$seurat_clusters)

combined_objects$K27ac_ac@meta.data$Annotation <- cluster_to_celltype[clusters] %>% sub("_1$|_2$", "", .)

# transfer annotation to other modalities
anno_K27ac <- data.frame(ID = Idents(combined_objects$K27ac_ac) %>% names(),
                         Annotation = combined_objects$K27ac_ac@meta.data$Annotation)

anno_ATAC <- data.frame(ID = Idents(combined_objects$ATAC_aa) %>% names())
anno_ATAC <- left_join(anno_ATAC, anno_K27ac)
combined_objects$ATAC_aa@meta.data$Annotation <- anno_ATAC$Annotation


anno_me <- data.frame(ID = Idents(combined_objects$K27me3_me) %>% names())
anno_me <- left_join(anno_me, anno_K27ac)
combined_objects$K27me3_me@meta.data$Annotation <- anno_me$Annotation

# export annotation and cellID for visualization of whole dataset
df_anno <-data.frame(CellID = paste0(REGION, rownames(combined_objects$K27ac_ac@meta.data)),
                     Annoation = combined_objects$K27ac_ac@meta.data$Annotation)

write.table(df_anno, paste0(RESULT_DIR_FINAL, "/nCT_", REGION, "_Annotation_using_peaks_after_reseq.txt"), sep = "\t", quote = F)


# ---- 14. run ChromVAR --------------------------------------------------------
pfm <- getMatrixSet(
  x = "JASPAR2024.sqlite3", # https://jaspar.elixir.no/downloads/
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# Add motifs
DefaultAssay(combined_objects$K27ac_ac) <- "peaks"

combined_objects$K27ac_ac <- AddMotifs(
  object = combined_objects$K27ac_ac,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

combined_objects$K27ac_ac <- RunChromVAR(
  object = combined_objects$K27ac_ac,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  new.assay.name = "chromvar"
)

saveRDS(combined_objects, paste0(RESULT_DIR_FINAL, "/", REGION, "_combined_SO_ChromVAR.rds"))
