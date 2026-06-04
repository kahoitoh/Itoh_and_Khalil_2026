#!/usr/bin/env Rscript

# ==============================================================================
# Classify AP-1-associated eRegulon enhancers by cell-type-specific ATAC activity
# ------------------------------------------------------------------------------
# Purpose:
#   This script quantifies ATAC accessibility at AP-1-associated eRegulon
#   enhancers identified by SCENIC+ and classifies them as cell-type-specific or
#   common enhancers.
#
# Main steps:
#   1. Load a processed Seurat object containing ATAC fragments and chromVAR scores.
#   2. Import AP-1-associated eRegulon regions extracted from SCENIC+ results.
#   3. Quantify ATAC fragment counts over AP-1 eRegulon enhancer regions.
#   4. Create a new chromatin assay for SCENIC+ enhancer regions and normalize it
#      using TF-IDF.
#   5. Select Cond4h cells and classify them as AP-1-high or AP-1-low based on
#      chromVAR AP-1 motif activity within each cell type.
#   6. Calculate mean enhancer accessibility for each cell type and AP-1 High/Low
#      group.
#   7. Annotate each enhancer by the SCENIC+ eRegulon cell type in which it was
#      detected.
#   8. Classify enhancers as cell-type-specific or common based on accessibility
#      margin and low-background accessibility thresholds.
#
# Inputs:
#   - SO_allCells_<REGION>.rds
#   - AP1_eRegulon_<REGION>.csv
#   - 10x Multiome ATAC fragment file
#
# Outputs:
#   - Classed_eRegulon_<REGION>.txt
#
# Notes:
#   - AP-1 motif ID: MA1141.2.
#   - AP-1-high cells are defined as the top 20% of chromVAR AP-1 scores within
#     each cell type.
#   - AP-1-low cells are defined as the bottom 20% of chromVAR AP-1 scores within
#     each cell type.
#   - For mPFC, excitatory subclasses are merged as Exc, and Sst/Pvalb cells are
#     merged as Inh.
#   - Enhancers are classified as Specific when AP-1-high accessibility is higher
#     than other cell types by at least delta_for_diff and low-state accessibility
#     is below the cell-type-specific 70th percentile.
# ==============================================================================

library(dplyr)
library(tibble)
library(Seurat)
library(Signac)
library(GenomicRanges)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(stringr)

# config -----------------------------------------------------------------------
REGION <- c("BLA", "Hippo", "mPFC")[3]
PATH_TO_SEURAT_OBJECT <- "path/to/your/SeuratObject"
PATH_TO_YOUR_FRAGMENT_FILE <- "path/to/your/atac_fragments.tsv.gz"
FINAL_RES_DIR <- "path/to/final/res"

FC_CONFIG <- if(REGION == "BLA"){
  tibble::tribble(
    ~celltype_out, ~celltype_col,
    "LA_Chst9",    "LA_Chst9",
    "BA",          "BA",
    "Sst",         "BLA-Sst"
  )
}else if(REGION == "Hippo"){
  tibble::tribble(
    ~celltype_out, ~celltype_col,
    "DG",    "DG",
    "CA1",          "CA1",
    "CA3",         "CA3"
  )
}else if(REGION == "mPFC"){
  tibble::tribble(
    ~celltype_out, ~celltype_col,
    "Exc",    "Exc",
    "Inh",         "Inh"
  )
}

FC_CONFIG <- FC_CONFIG %>%
  mutate(
    high_col = paste0(celltype_col, "-High"),
    low_col  = paste0(celltype_col, "-Low")
  )

# import data ------------------------------------------------------------------
SO <- readRDS(paste0(PATH_TO_SEURAT_OBJECT, "/SO_allCells_", REGION, ".rds"))

# Calculate ATAC signals within enhancers identified in Scenicplus analysis
Scenicplus_enhancers <- read.csv(paste0(FINAL_RES_DIR, "/AP1_eRegulon_", REGION, ".csv"))
Scenicplus_enhancers <- Scenicplus_enhancers$Region %>% unique()
Scenicplus_enhancers_bed <- data.frame(chr = sub(":.+", "", Scenicplus_enhancers),
                                       start = sub("-[0-9]+", "", Scenicplus_enhancers) %>% sub(".+:", "", .) %>% as.numeric(),
                                       end = sub(".+-", "", Scenicplus_enhancers) %>% as.numeric())

peaks_gr <- GRanges(
  seqnames = Scenicplus_enhancers_bed$chr,
  ranges   = IRanges(start = Scenicplus_enhancers_bed$start + 1,
                     end   = Scenicplus_enhancers_bed$end)
)

enhancer_counts <- FeatureMatrix(
  fragments = Fragments(SO),
  features = peaks_gr,
  cells = colnames(SO)
)

annotation <- GetGRangesFromEnsDb(EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"

SO[["Scenicplus_Enhancers"]] <- CreateChromatinAssay(
  counts = enhancer_counts,
  fragments = PATH_TO_YOUR_FRAGMENT_FILE,
  annotation = annotation
)

DefaultAssay(SO) <- "Scenicplus_Enhancers"
SO <- RunTFIDF(SO)

# perform Trapping -------------------------------------------------------------
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

# Calculate ATAC signal in AP-1 eRegulon enhancers in TRAP High and Low cells --
SO_Cond_High_Low@meta.data$Cell_Cond <- paste0(SO_Cond_High_Low@meta.data$Annotation, "-", SO_Cond_High_Low@meta.data$High_low)

if(REGION == "mPFC"){
  SO_Cond_High_Low@meta.data$Cell_Cond <- sub("L2_3_IT|L4_5_IT|L5_ET|L5_NP|L6_CT|L6_IT", "Exc", SO_Cond_High_Low@meta.data$Cell_Cond) %>% sub("Sst|Pvalb", "Inh", .)
}

Scenicplus_Enhancers_data <- GetAssayData(SO_Cond_High_Low, assay = "Scenicplus_Enhancers", layer = "data")

cell_cond <- SO_Cond_High_Low@meta.data[colnames(Scenicplus_Enhancers_data), "Cell_Cond"]
count_mean_by_cond <- sapply(
  split(seq_along(cell_cond), cell_cond),
  function(idx) {
    Matrix::rowMeans(Scenicplus_Enhancers_data[, idx, drop = FALSE])
  }
) %>% as.data.frame() 

# take log2 FC of ATAC-seq levels ----------------------------------------------
count_mean_selected <- count_mean_by_cond[,c(FC_CONFIG$low_col, FC_CONFIG$high_col)]
count_mean_selected$PeakName <- rownames(count_mean_selected)

# Annotate peaks ---------------------------------------------------------------
Scenicplus_enhancers <- read.csv(paste0(FINAL_RES_DIR, "/AP1_eRegulon_", REGION, ".csv"))

Scenicplus_enhancers$chr <- sub(":.+", "", Scenicplus_enhancers$Region)
Scenicplus_enhancers$start <- sub("-[0-9]+$", "", Scenicplus_enhancers$Region) %>% sub(".+:", "", .) %>% as.numeric()
Scenicplus_enhancers$end <- sub(".+-", "", Scenicplus_enhancers$Region)

Scenicplus_enhancers$PeakName <- paste0(Scenicplus_enhancers$chr, "-", Scenicplus_enhancers$start + 1, "-", Scenicplus_enhancers$end)

Scenicplus_enhancers_summary <- Scenicplus_enhancers[,c("PeakName", "CellType")] %>% table()

anno <- c()
for (i in seq_along(Scenicplus_enhancers_summary[,1])){
  
  anno <- c(anno, which(Scenicplus_enhancers_summary[i,] != 0) %>% names() %>% str_c(., collapse = "-"))
  
}

df_anno <- data.frame(PeakName = Scenicplus_enhancers_summary %>% rownames(),
                      Annotation = anno)
df_anno <- df_anno[order(df_anno$Annotation),]

count_mean_selected <- left_join(df_anno, count_mean_selected)

# Calculate log2 FC of ATAC signal in AP-1 enhancers among cell types ----------
delta_for_diff <- 0.1

if(REGION %in% c("BLA", "Hippo")){
  
  count_mean_selected <- count_mean_selected %>%
    rowwise() %>%
    mutate(
      margin = case_when(
        Annotation == FC_CONFIG$celltype_out[1]  ~ .data[[FC_CONFIG$high_col[1]]]  - max(c(.data[[FC_CONFIG$high_col[2]]], .data[[FC_CONFIG$high_col[3]]]), na.rm = TRUE),
        Annotation == FC_CONFIG$celltype_out[2]  ~ .data[[FC_CONFIG$high_col[2]]]  - max(c(.data[[FC_CONFIG$high_col[1]]], .data[[FC_CONFIG$high_col[3]]]), na.rm = TRUE),
        Annotation == FC_CONFIG$celltype_out[3]  ~ .data[[FC_CONFIG$high_col[3]]]  - max(c(.data[[FC_CONFIG$high_col[1]]], .data[[FC_CONFIG$high_col[2]]]), na.rm = TRUE),
        TRUE ~ NA_real_
      ),
      new_val = case_when(
        Annotation == FC_CONFIG$celltype_out[1] ~ .data[[FC_CONFIG$low_col[1]]],
        Annotation == FC_CONFIG$celltype_out[2] ~ .data[[FC_CONFIG$low_col[2]]],
        Annotation == FC_CONFIG$celltype_out[3] ~ .data[[FC_CONFIG$low_col[3]]],
        TRUE ~ NA_real_
      ),
      pass_thresh = margin >= delta_for_diff
    ) %>%
    ungroup() %>% 
    group_by(Annotation) %>%
    mutate(
      q70_new = quantile(new_val, probs = 0.7, na.rm = TRUE),
      pass_q70 = new_val < q70_new,
      pass = pass_thresh & pass_q70
    ) %>%
    ungroup()
  
}else if(REGION == "mPFC"){
  
  count_mean_selected <- count_mean_selected %>%
    rowwise() %>%
    mutate(
      margin = case_when(
        Annotation == FC_CONFIG$celltype_out[1]  ~ .data[[FC_CONFIG$high_col[1]]]  - .data[[FC_CONFIG$high_col[2]]],
        Annotation == FC_CONFIG$celltype_out[2]  ~ .data[[FC_CONFIG$high_col[2]]]  - .data[[FC_CONFIG$high_col[1]]],
        TRUE ~ NA_real_
      ),
      new_val = case_when(
        Annotation == FC_CONFIG$celltype_out[1] ~ .data[[FC_CONFIG$low_col[1]]],
        Annotation == FC_CONFIG$celltype_out[2] ~ .data[[FC_CONFIG$low_col[2]]],
        TRUE ~ NA_real_
      ),
      pass_thresh = margin >= delta_for_diff
    ) %>%
    ungroup() %>% 
    group_by(Annotation) %>%
    mutate(
      q70_new = quantile(new_val, probs = 0.7, na.rm = TRUE),
      pass_q70 = new_val < q70_new,
      pass = pass_thresh & pass_q70
    ) %>%
    ungroup()
  
}

count_mean_selected <- count_mean_selected[which(count_mean_selected$pass_q70 == T),]
count_mean_selected$Enhancer_class <- case_when(count_mean_selected$pass == T~"Specific",
                                                count_mean_selected$pass_thresh == F~"Common")

write.table(count_mean_selected, 
            paste0(FINAL_RES_DIR, "/Classed_eRegulon_", REGION, ".txt"), quote = F, sep = "\t")

# Make bedfile of cell type-specific enhancers ---------------------------------
Sepecific_regulon <- count_mean_selected[which(count_mean_selected$Enhancer_class == "Specific"),]

for (ct in FC_CONFIG$celltype_out){
  
  df_bed <- data.frame(Chr = sub("-.+", "", Sepecific_regulon$PeakName[which(Sepecific_regulon$Annotation == ct)]) %>% sub("chr", "", .),
                       Start = sub("-[0-9]+", "", Sepecific_regulon$PeakName[which(Sepecific_regulon$Annotation == ct)]) %>% sub(".+-", "", .),
                       End = sub(".+-", "", Sepecific_regulon$PeakName[which(Sepecific_regulon$Annotation == ct)]))
  
  write.table(df_bed, 
              paste0(FINAL_RES_DIR, "/bed_SpecificEnhancer_", REGION, "_", ct, ".bed"), quote = F, sep = "\t", col.names = F, row.names = F)
  
}
