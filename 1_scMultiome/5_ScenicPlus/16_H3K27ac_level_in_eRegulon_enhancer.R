#!/usr/bin/env Rscript

# ==============================================================================
# H3K27ac signal at cell-type-specific AP-1 eRegulon enhancers
# ------------------------------------------------------------------------------
# This script projects SCENIC+ AP-1 eRegulon enhancer regions onto nanoCT H3K27ac
# data and compares enhancer-associated H3K27ac signal between AP-1/FOS
# chromVAR-high and chromVAR-low cell populations.
#
# The workflow is organized into two analyses:
#
#   1. Heatmap of enhancer-associated H3K27ac signal
#      - selects AP-1 eRegulon enhancers annotated as cell-type specific
#      - calculates the mean TF-IDF-normalized H3K27ac signal for each enhancer
#        class across AP-1/FOS-high and AP-1/FOS-low cells
#      - z-scores the mean signal within each enhancer class to highlight
#        relative enrichment patterns across cell populations
#
#   2. Permutation test for enhancer detection enrichment
#      - aggregates raw H3K27ac counts by cell-type / AP-1-FOS state / replicate
#      - tests whether each cell-type-specific enhancer class is detected more
#        frequently in the matching AP-1/FOS-high population
#      - uses random enhancer sets matched by overall detection frequency as the
#        empirical null distribution
#
# Inputs:
#   - <REGION>_combined_SO_ChromVAR_ATAC.rds
#   - AP1_eRegulon_<REGION>.csv
#   - Classed_eRegulon_<REGION>.txt
#
# Outputs:
#   - Heatmap_eRegulon_H3K27ac_Scale_within_Enhancer_<REGION>.png
#   - Permutation_res_list and BH-adjusted empirical p-values
#
# Notes:
#   - Set REGION to BLA, Hippo, or mPFC before running the script.
#   - For mPFC, neuronal annotations are collapsed into Exc and Inh groups.
#   - High/Low groups are defined from the top/bottom 20% of AP-1/FOS chromVAR
#     scores within each annotated cell type, followed by selection of Cond cells.
#   - The heatmap summarizes relative TF-IDF-normalized H3K27ac signal.
# ==============================================================================

library(Seurat)
library(Signac) 
library(EnsDb.Mmusculus.v79)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(stringr)
library(purrr)

# config -----------------------------------------------------------------------
PATH_TO_eREGULON <- "D:/in_vivo_Multiome_EngramProject/Github/1_scMultiome/5_ScenicPlus/data"
FINAL_RES_PATH <- "D:/in_vivo_Multiome_EngramProject/Github/2_multi-nanoCT/1_Preprocess/data"
REGION <- c("BLA", "Hippo", "mPFC")[1]
RESULT_DIR_FINAL <- paste0("D:/in_vivo_nanoCT_2502/nCT_", REGION, "_baseres_rm_HC1_ATAC_Cells_Having_All_modalities")

CELLTYPE_to_CHOOSE <- if(REGION == "BLA"){
  
  list("LA", "BA",  "Sst")
  
}else if(REGION == "Hippo"){
  
  list("DG", "CA1", "CA3")
  
}else if(REGION == "mPFC"){
  
  list("Exc", "Inh")
  
}

# Helper -----------------------------------------------------------------------
get_ac_levels_CellTypeSpecificEnhancers <- function(ct){
  
  # take cell type specific enhancers
  Specific_eRegulon_ct <- Specific_eRegulon$PeakName[grep(ct, Specific_eRegulon$Annotation)]
  index_enhancer <- which(rownames(Scenicplus_Enhancers_data) %in% Specific_eRegulon_ct)
  
  # calculate ac levels in High and Low populations in different cell types
  df_mean_list <- list()
  for (ct2 in CELLTYPE_to_CHOOSE){
    
    index_celltype <- which(SO_TRAPed_ac@meta.data$Annotation == ct2)
    
    Scenicplus_Enhancers_data_selected <- Scenicplus_Enhancers_data[index_enhancer,index_celltype]
    
    Scenicplus_Enhancers_mean <- Scenicplus_Enhancers_data_selected %>% colMeans() %>% as.data.frame()
    colnames(Scenicplus_Enhancers_mean)[1] <- "SignalInt"
    Scenicplus_Enhancers_mean$high_low <- SO_TRAPed_ac@meta.data$High_low[which(SO_TRAPed_ac@meta.data$Annotation == ct2)]
    
    df_mean_list[[ct2]] <- Scenicplus_Enhancers_mean %>%
      group_by(high_low) %>%
      summarise(
        mean = mean(SignalInt, na.rm = TRUE),
        n_genes = n(),
        .groups = "drop"
      ) %>% 
      mutate(CellType = ct2)
  }
  
  res <- do.call(rbind, df_mean_list)
  res$Cell_Cond <- paste0(res$CellType, "-", res$high_low)
  res$Enhancers <- ct
  
  res
  
}

perm_test_enhancer <- function(target_class,
                               target_celltype,
                               
                               counts,
                               peak_class,
                               sample_celltype,
                               
                               stat = c("diff_detect", "logOR"),
                               nperm = 5000,
                               seed = 1,
                               stratify_by = c("none", "overall_detect_bin"),
                               nbins = 10) {
  
  stat <- match.arg(stat)
  stratify_by <- match.arg(stratify_by)
  
  stopifnot(ncol(counts) == length(sample_celltype))
  stopifnot(nrow(counts) == length(peak_class))
  
  set.seed(seed)
  
  # detection matrix (peak x sample): 1 if detected
  det <- (counts > 0) * 1L
  
  # target samples vs others
  targ <- sample_celltype == target_celltype
  oth  <- !targ
  if (sum(targ) < 1 || sum(oth) < 1) stop("Need >=1 target and >=1 other sample")
  
  # observed set index
  idx_set <- peak_class == target_class
  m <- sum(idx_set)
  if (m < 2) stop("Set size < 2; cannot test reliably")
  
  # statistic calculator for a given boolean peak-set index
  calc_stat <- function(idx) {
    # aggregate across peaks and samples (set-level)
    a <- sum(det[idx, targ, drop=FALSE])  # detected in target
    b <- sum(1L - det[idx, targ, drop=FALSE])  # not detected in target
    c <- sum(det[idx, oth,  drop=FALSE])  # detected in others
    d <- sum(1L - det[idx, oth,  drop=FALSE])
    
    if (stat == "diff_detect") {
      # detection rate difference (target - other)
      p1 <- a / (a + b)
      p2 <- c / (c + d)
      return(p1 - p2)
    } else {
      # log odds ratio with +0.5 pseudocount (Haldane-Anscombe) for stability
      return(log(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))))
    }
  }
  
  obs <- calc_stat(idx_set)
  
  # permutation scheme
  if (stratify_by == "none") {
    # shuffle class labels across all peaks
    perm_stats <- replicate(nperm, {
      perm_idx <- sample.int(nrow(counts), m, replace = FALSE)
      idx <- rep(FALSE, nrow(counts))
      idx[perm_idx] <- TRUE
      calc_stat(idx)
    })
  } else {
    # stratified shuffle by overall detection frequency (controls "easy-to-detect" peaks)
    overall <- rowMeans(det)  # per-peak detection frequency across all samples
    # bins <- cut(overall, breaks = quantile(overall, probs = seq(0, 1, length.out = nbins + 1), na.rm=TRUE),
    #             include.lowest = TRUE)
    
    qs <- quantile(
      overall,
      probs = seq(0, 1, length.out = nbins + 1),
      na.rm = TRUE,
      names = FALSE
    )
    
    qs <- unique(qs)
    
    if (length(qs) < 2) {
      bins <- factor(rep("all", length(overall)))
    } else {
      bins <- cut(
        overall,
        breaks = qs,
        include.lowest = TRUE
      )
    }
    
    
    # how many set peaks in each bin?
    tab_need <- table(bins[idx_set])
    
    # list peak indices by bin
    peaks_by_bin <- split(seq_len(nrow(counts)), bins)
    
    perm_stats <- replicate(nperm, {
      chosen <- integer(0)
      for (bn in names(tab_need)) {
        k <- as.integer(tab_need[[bn]])
        pool <- peaks_by_bin[[bn]]
        if (length(pool) < k) stop("Bin pool too small; reduce nbins or use stratify_by='none'")
        chosen <- c(chosen, sample(pool, k, replace = FALSE))
      }
      idx <- rep(FALSE, nrow(counts))
      idx[chosen] <- TRUE
      calc_stat(idx)
    })
  }
  
  # empirical p (one-sided: enriched in target => larger stat)
  p_emp <- (1 + sum(perm_stats >= obs)) / (nperm + 1)
  
  list(
    target_class = target_class,
    stat = stat,
    stratify_by = stratify_by,
    observed = obs,
    perm_mean = mean(perm_stats),
    perm_sd = sd(perm_stats),
    p_empirical = p_emp,
    nperm = nperm,
    set_size = m
  )
}



# Load data --------------------------------------------------------------------
# import seurat object of multi-nanoCT
combined_objects <- readRDS(paste0(RESULT_DIR_FINAL, "/", REGION, "_combined_SO_ChromVAR_ATAC.rds"))
if(REGION == "mPFC"){
  combined_objects$ATAC_aa@meta.data$Annotation <- sub("L2_3_IT|L4_5_IT|L5_ET|L5_NP|L6_CT|L6_IT", "Exc", combined_objects$ATAC_aa@meta.data$Annotation)
  combined_objects$ATAC_aa@meta.data$Annotation <- sub("Sst|PV", "Inh", combined_objects$ATAC_aa@meta.data$Annotation)
}


# import full eRegulon
eRegulon_full <- read.csv(paste0(PATH_TO_eREGULON, "/AP1_eRegulon_", REGION, ".csv"))
Scenicplus_enhancers <- eRegulon_full$Region %>% unique()
Scenicplus_enhancers_bed <- data.frame(chr = sub(":.+", "", Scenicplus_enhancers),
                                       start = sub("-[0-9]+", "", Scenicplus_enhancers) %>% sub(".+:", "", .) %>% as.numeric(),
                                       end = sub(".+-", "", Scenicplus_enhancers) %>% as.numeric())
peaks_gr <- GRanges(
  seqnames = Scenicplus_enhancers_bed$chr,
  ranges   = IRanges(start = Scenicplus_enhancers_bed$start + 1,
                     end   = Scenicplus_enhancers_bed$end)
)

# import specific eRegulon enhancers
Specific_eRegulon <- read.csv(paste0(PATH_TO_eREGULON, "/Classed_eRegulon_", REGION, ".txt"),
                                     header = T, sep = "\t")
Specific_eRegulon <- Specific_eRegulon[which(Specific_eRegulon$Enhancer_class == "Specific"),]

# import annotation 
annotation <- GetGRangesFromEnsDb(EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"

# Calculate H3K27ac levels in eRegulon enhancers -------------------------------
enhancer_counts <- FeatureMatrix(
  fragments = Fragments(combined_objects$K27ac_ac),
  features = peaks_gr,
  cells = colnames(combined_objects$K27ac_ac)
)

combined_objects$K27ac_ac[["Scenicplus_Enhancers"]] <- CreateChromatinAssay(
  counts = enhancer_counts,
  fragments = Fragments(combined_objects$K27ac_ac),
  annotation = annotation
)

DefaultAssay(combined_objects$K27ac_ac) <- "Scenicplus_Enhancers"
combined_objects$K27ac_ac <- RunTFIDF(combined_objects$K27ac_ac)

combined_objects$K27ac_ac@meta.data$FOS_MotifScore <- combined_objects$K27ac_ac@assays$chromvar@data["MA1141.2",]

# perform Trapping -------------------------------------------------------------
chromvar_count_selected <- combined_objects$ATAC_aa@assays$chromvar$data[which(rownames(combined_objects$ATAC_aa@assays$chromvar$data) == "MA1141.2"),]
df_chromvar_count <- data.frame(ChromVAR_count = chromvar_count_selected)
df_chromvar_count$CellID <- rownames(df_chromvar_count)
df_chromvar_count$Celltype <- combined_objects$ATAC_aa@meta.data$Annotation
df_chromvar_count$Sample <- combined_objects$ATAC_aa@meta.data$sample %>% sub("1|2", "", .)

df_chromvar_count_High <- df_chromvar_count %>% 
  group_by(., Celltype) %>% 
  mutate(q80 = quantile(ChromVAR_count, 0.8, na.rm = TRUE)) %>%
  dplyr::filter(ChromVAR_count >= q80 &
                  Sample == "Cond")

df_chromvar_count_Low <- df_chromvar_count %>% 
  group_by(., Celltype) %>% 
  mutate(q20 = quantile(ChromVAR_count, 0.2, na.rm = TRUE)) %>%
  dplyr::filter(ChromVAR_count <= q20 &
                  Sample == "Cond")

combined_objects$K27ac_ac@meta.data$CellID <- rownames(combined_objects$K27ac_ac@meta.data)

SO_TRAPed_ac <- subset(combined_objects$K27ac_ac, CellID %in% c(df_chromvar_count_High$CellID, df_chromvar_count_Low$CellID))
SO_TRAPed_ac@meta.data$High_low <- case_when(SO_TRAPed_ac@meta.data$CellID %in% df_chromvar_count_Low$CellID~"Low",
                                                    TRUE~"High")
SO_TRAPed_ac@meta.data$Cell_Cond <- paste0(SO_TRAPed_ac@meta.data$Annotation, "-", SO_TRAPed_ac@meta.data$High_low)

if(REGION == "mPFC"){
  SO_TRAPed_ac@meta.data$Annotation <- sub("L2_3_IT|L4_5_IT|L5_ET|L5_NP|L6_CT|L6_IT", "Exc", SO_TRAPed_ac@meta.data$Annotation)
  SO_TRAPed_ac@meta.data$Annotation <- sub("Sst|PV", "Inh", SO_TRAPed_ac@meta.data$Annotation)
}


# get H3K27ac levels in eRegulon enhancers in TRAPed neurons -------------------
Scenicplus_Enhancers_data <- GetAssayData(SO_TRAPed_ac, assay = "Scenicplus_Enhancers", layer = "data")

# filter the result for cell type-specific eRegulon in different cell types
# and take mean ----------------------------------------------------------------
Enhancer_ac_levels_list <- lapply(CELLTYPE_to_CHOOSE, get_ac_levels_CellTypeSpecificEnhancers)
Enhancer_ac_levels <- do.call(rbind, Enhancer_ac_levels_list)

# normalize among enhanceres
Enhancer_ac_levels <- Enhancer_ac_levels %>%
  group_by(Enhancers) %>%
  mutate(mean_z = as.numeric(scale(mean)),
         Enhancers_Cond = paste0(Enhancers, "-", high_low)) %>%  # グループ内Z
  ungroup()

Enhancer_ac_levels$Enhancers <- factor(Enhancer_ac_levels$Enhancers, CELLTYPE_to_CHOOSE[length(CELLTYPE_to_CHOOSE):1])

g_Specific_Z_scale_within_Enhancer <- ggplot(Enhancer_ac_levels, aes(x =Cell_Cond, y = Enhancers, fill = mean_z))+
  theme_minimal()  +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_fill_viridis_c(option = "plasma") +
  ylab("Enhancer") +
  xlab("Cells")

ggsave(paste0(FINAL_RES_PATH, "/Heatmap_eRegulon_H3K27ac_Scale_within_Enhancer_", REGION, ".png"), g_Specific_Z_scale_within_Enhancer,
       width = 5, height = 3)

# perform Permutation test for H3K27ac levels in AP-1 eRegulon enhancer --------

# calculate sum ac levels in enhancers per cell type
cell_cond <- paste0(SO_TRAPed_ac@meta.data$Cell_Cond, "-", sub("Cond", "", SO_TRAPed_ac@meta.data$orig.ident))

if(REGION == "mPFC"){
  cell_cond <- sub("L2_3_IT|L4_5_IT|L5_ET|L5_NP|L6_CT|L6_IT", "Exc", cell_cond)
  cell_cond <- sub("Sst|PV|Sst_PV", "Inh", cell_cond)
}

Scenicplus_Enhancers_count <- GetAssayData(SO_TRAPed_ac, assay = "Scenicplus_Enhancers", layer = "counts")

Enhancer_ac_Sum_by_Cell_cond <- sapply(
  split(seq_along(cell_cond), cell_cond),
  function(idx) {
    Matrix::rowSums(Scenicplus_Enhancers_count[, idx, drop = FALSE])
  }
) %>% as.data.frame() 

Enhancer_ac_Sum_by_Cell_cond <- Enhancer_ac_Sum_by_Cell_cond[,grep(paste(CELLTYPE_to_CHOOSE, collapse = "|"), colnames(Enhancer_ac_Sum_by_Cell_cond))]

counts <- as.matrix(Enhancer_ac_Sum_by_Cell_cond)

sample_celltype <- colnames(counts) %>% sub("-1|-2", "", .)

# get ac levels in specific enhancers
index_peakname <- which(rownames(counts) %in% Specific_eRegulon$PeakName)
counts_selected <- counts[index_peakname,]

tmp_counts_peakname <- data.frame(PeakName = rownames(counts_selected))
peak_anno_aligned <- left_join(tmp_counts_peakname,
                               Specific_eRegulon)

peak_class <- peak_anno_aligned$Annotation %>% sub("_Chst9", "", .)
names(peak_class) <- peak_anno_aligned$PeakName

# perform permutation test
Permutation_res_list <- map2(CELLTYPE_to_CHOOSE,
                             paste0(CELLTYPE_to_CHOOSE, "-High"),
                             perm_test_enhancer,
                             counts_selected,
                             peak_class,
                             sample_celltype,
                             "logOR",
                             5000,
                             3,
                             "overall_detect_bin",
                             5)

lapply(Permutation_res_list, function(res){res$p_empirical})

# Adjust p values
p_emp_adj <- p.adjust(lapply(Permutation_res_list, function(res){res$p_empirical}) %>% unlist(), method = "BH")
names(p_emp_adj) <- CELLTYPE_to_CHOOSE
p_emp_adj 
