#!/usr/bin/env Rscript

# ==============================================================================
# DEG identification for ChromTRAP benchmark
# a. ChromTRAP Cond High versus ChromTRAP HC Low DEG 
# ------------------------------------------------------------------------------
# This script identifies DEGs between ChromTRAP-defined Cond High and HC Low
# excitatory neurons. HC Low cells are defined as neurons below the 10th
# percentile of chromVAR activity for the selected motif in HC samples.
# Cond High cells are defined as neurons above the 90th percentile in the
# Cond 4 h sample. Pseudobulk RNA counts are generated with Seurat
# AggregateExpression and DEGs are identified using DESeq2.
#
# Edit the CONFIG section before running.
# ==============================================================================

library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(Seurat)
library(DESeq2)

# Config------------------------------------------------------------------------
region <- "BLA" # or "Hippo"

SO_VGlut <- readRDS("D:/in_vivo_Multiome_EngramProject/BaseAnalysis_BLA/6_ATAC_ChromVAR/SO_VGlut_chromVAR.rds") # seurat object containing Exc neurons in whole nuclei sequencing data

CellTypes_selected <- if(region == "Hippo"){
  c("DG_Neu_1", "DG_Neu_2", "CA1_Neu", "CA3_Neu", "CA3_Neu_St18")
  }else if(region == "BLA"){
    c("LA_Neu_Rorb", "LA_Neu_Chst9", "BA_Neu")
  }

finalres_dir <- "path_to_final_res"

# ChromTRAP --------------------------------------------------------------------
SO_VGlut <- subset(SO_VGlut, SCT_snn_res.0.15_Anno %in% CellTypes_selected)

# get neurons  below the 10th percentile in the HC: "HC Low"
SO_HC <- subset(SO_VGlut, Sample %in% c("HC_1", "HC_2"))

chromvar_count_selected <- SO_HC@assays$chromvar$data[which(rownames(SO_HC@assays$chromvar$data) == "MA1141.2"),]
df_chromvar_count <- data.frame(ChromVAR_count = chromvar_count_selected)

df_chromvar_count$CellID <- rownames(df_chromvar_count)

df_chromvar_count_q10 <- df_chromvar_count %>% 
  mutate(q10 = quantile(ChromVAR_count, 0.1, na.rm = TRUE)) %>%
  dplyr::filter(ChromVAR_count <= q10)

# get neurons above the 90th percentile in the Cond 4 h group: "Cond High"
SO_Cond4h <- subset(SO_VGlut, Sample %in% c("Cond4h_1")) # use only rep1 to align with FOS-FACS data

chromvar_count_selected <- SO_Cond4h@assays$chromvar$data[which(rownames(SO_Cond4h@assays$chromvar$data) == "MA1141.2"),]
df_chromvar_count <- data.frame(ChromVAR_count = chromvar_count_selected)

df_chromvar_count$CellID <- rownames(df_chromvar_count)

df_chromvar_count_q90 <- df_chromvar_count %>% 
  mutate(q90 = quantile(ChromVAR_count, 0.9, na.rm = TRUE)) %>%
  dplyr::filter(ChromVAR_count >= q90)

# filter seurat object for "HC Low" and "Cond High"
SO_VGlut@meta.data$CellID <- rownames(SO_VGlut@meta.data)
SO_TRAP <- subset(SO_VGlut, CellID %in% c(df_chromvar_count_q90$CellID, df_chromvar_count_q10$CellID))


# Take pseuodbulk for each condition -------------------------------------------
pb <- AggregateExpression(SO_TRAP, group.by = c("Sample"), assays = "RNA",
                          slot = "counts")$RNA

samples <- colnames(pb)

coldata <- data.frame(Condition = sub("-1|-2", "", samples),
                      Celltype_condition = samples)

# Perform DEseq2 ---------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = round(pb),
                              colData = coldata,
                              design = ~ Condition) 

dds <- DESeq(dds)  # size factor estimation + normalization

norm_counts <- counts(dds, normalized = TRUE)  

res <- results(dds, contrast=c("Condition","Cond4h","HC"))

res_df <- res %>% as.data.frame()

res_df$Class_label <- case_when(res_df$log2FoldChange > 1 & res_df$padj < 0.07~rownames(res_df),
                                TRUE~"")

res_df$Class_label_selected <- case_when(res_df$Class_label %in% c("Elmo1", "Vgf", "Bdnf", "Sorcs3", "Sorcs1","Nptx2", "Scg2")~res_df$Class_label,
                                         TRUE~"")

res_df$Class <- case_when(res_df$Class_label != ""~"Upregulated",
                          TRUE~"other")

g_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Class, label = Class_label_selected, alpha = Class)) + 
  geom_point(alpha = 0.2) +
  geom_text_repel(max.overlaps = Inf) +
  theme_cowplot() +
  geom_vline(xintercept = 1, linetype = 3) +
  geom_hline(yintercept = -log10(0.07), linetype = 3) +
  scale_color_manual(values = c("other" = "gray",
                                "Upregulated" = "deeppink")) +
  scale_alpha_manual(values = c("other" = 0.1,
                                "Upregulated" = 0.5)) +
  labs(subtitle = paste0("Cond High: ", length(which(SO_TRAP@meta.data$Sample == "Cond4h_1")), " cells\n",
                         "HC Low rep 1: ", length(which(SO_TRAP@meta.data$Sample == "HC_1")), " cells\n",
                         "HC Low rep 2: ", length(which(SO_TRAP@meta.data$Sample == "HC_2")), " cells"))

ggsave(paste0(finalres_dir, "/", region, "_Volcano_CondHigh_vs_HCLow.png"), g_volcano, width = 5, height = 4)
ggsave(paste0(finalres_dir, "/", region, "_Volcano_CondHigh_vs_HCLow.pdf"), g_volcano, width = 5, height = 4)

write.table(res_df, paste0(finalres_dir, "/", region, "_DEseqres_CondHigh_vs_HCLow.txt"),
            quote = F, sep = "\t")

