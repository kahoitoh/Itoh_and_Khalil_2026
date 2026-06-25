#!/usr/bin/env Rscript

# ==============================================================================
# DEG identification for ChromTRAP benchmark
# b. FOS-FACS versus ChromTRAP HC Low DEG analysis
# ------------------------------------------------------------------------------
# This script identifies DEGs between FOS-FACS-sorted Cond 4 h neurons and
# ChromTRAP-defined HC Low neurons, then compares the resulting DEG set with
# DEGs from the Cond High versus HC Low ChromTRAP benchmark.
#
# HC Low cells are defined as HC neurons below the 10th percentile of chromVAR
# activity for the selected motif. Pseudobulk RNA counts are generated with
# Seurat AggregateExpression, mitochondrial/ribosomal/background genes are
# removed, and DEGs are identified with DESeq2.
#
# Edit the CONFIG section before running. 
# ==============================================================================

library(Seurat)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(eulerr)

# Config------------------------------------------------------------------------
region <- "BLA" # or "Hippo"

PATH_TO_SO_FACS_MO <- "path/to/your/SO_FACS_MO"
SO <- readRDS(paste0(finalres_dir, "/", region, "_", "SO_chromVAR.rds")) # from 1_scMultiome/1_Preprocess/3_PreProcess_FOS_FACS_scMultiome.R

CellTypes_selected <- if(region == "Hippo"){
  c("DG_1", "DG_2", "CA1", "CA3_Col6a1", "CA3_Kcnq5")
}else if(region == "BLA"){
  c("BA", "LA")
}

PATH_TO_DEseq2res_FACS_MO <- "path/to/your/SO_FACS_MO/DESeq2/res"
DEseq2res_CondHigh_vs_HCLow <- read.table(paste0(PATH_TO_DEseq2res_FACS_MO, "/", region, "_DEseqres_CondHigh_vs_HCLow.txt"), 
                                          header = T, sep = "\t") # output from 1_scMultiome/4_DEG_to_benchmark_ChromTRAP/1_DEG_to_benchmark_ChromTRAP_CondHigh_vs_HCLow.R

finalres_dir <- "path/to/final/res" 

# chromTRAP---------------------------------------------------------------------
SO_ExcNeu <- subset(SO, Annotation %in% CellTypes_selected)

SO_HC <- subset(SO_ExcNeu, Sample %in% c("HC_OptiPrep_1", "HC_OptiPrep_2"))

chromvar_count_selected <- SO_HC@assays$chromvar$data[which(rownames(SO_HC@assays$chromvar$data) == "MA1141.2"),]
df_chromvar_count <- data.frame(ChromVAR_count = chromvar_count_selected)
df_chromvar_count$CellID <- rownames(df_chromvar_count)

df_chromvar_count_q10 <- df_chromvar_count %>% 
  mutate(q10 = quantile(ChromVAR_count, 0.1, na.rm = TRUE)) %>%
  dplyr::filter(ChromVAR_count <= q10)

SO_ExcNeu@meta.data$CellID <- rownames(SO_ExcNeu@meta.data)

# filter seurat object
SO_FACS_HCLow <- subset(SO_ExcNeu, Sample == "Cond4h_FACS" | CellID %in% df_chromvar_count_q10$CellID)

# perform DEseq2 ---------------------------------------------------------------
pb <- AggregateExpression(SO_FACS_HCLow, group.by = c("Sample"), assays = "RNA",
                          slot = "counts")$RNA

# remove mitochondrial and ribosomal genes
is_mito <- rownames(pb)[grep("^Mt-", rownames(pb), ignore.case=TRUE)]
is_ribo <- rownames(pb)[grep("^Rpl|^Rps", rownames(pb))]
blocklist <- c(is_mito, is_ribo,
               "Malat1","Neat1","Ubb","Ubc","Hsp90aa1","Hsp90ab1","Hspa1a","Hspa1b")

pb_rm_block <- pb[-which(rownames(pb) %in% blocklist),]

samples <- colnames(pb_rm_block)

coldata <- data.frame(Celltype = sub("-Cond4h-FACS|-HC-OptiPrep", "", colnames(pb_rm_block)),
                      Condition = sub("-1|-2", "", colnames(pb_rm_block)) %>% sub(".+-", "", .),
                      Celltype_condition = colnames(pb_rm_block) %>% as.factor())

dds <- DESeqDataSetFromMatrix(countData = round(pb_rm_block),
                              colData = coldata,
                              design = ~ Condition) 
dds <- DESeq(dds)  # size factor estimation + normalization
norm_counts <- counts(dds, normalized = TRUE)  

res <- results(dds, contrast=c("Condition","FACS","OptiPrep"))

res_df <- res %>% as.data.frame()

# import DEseq2 result of Cond4h High vs HC Low --------------------------------
DEseq2res_CondHigh_vs_HCLow_DEG <- rownames(DEseq2res_CondHigh_vs_HCLow)[which(DEseq2res_CondHigh_vs_HCLow$Class == "Upregulated")]


res_df$ARG_ChromTRAP_label <- case_when(rownames(res_df) %in% DEseq2res_CondHigh_vs_HCLow_DEG~rownames(res_df),
                                        TRUE~"")

res_df$ARG_ChromTRAP <- case_when(res_df$ARG_ChromTRAP_label != ""~T,
                                  TRUE~F)

res_df$ARG_ChromTRAP_FACS <- case_when(res_df$log2FoldChange > 1 & res_df$padj < 0.07 & res_df$ARG_ChromTRAP != ""~T,
                                       TRUE~F)

res_df$Class <- case_when(res_df$log2FoldChange > 1 & res_df$padj < 0.07~"Upregulated",
                          TRUE~"others")

res_df$Anno_threeColor <- case_when(res_df$ARG_ChromTRAP == T & res_df$Class == "Upregulated"~"DEG_both",
                                    res_df$ARG_ChromTRAP == F & res_df$Class == "Upregulated"~"DEG_FACS_only",
                                    res_df$ARG_ChromTRAP == T & res_df$Class == "others"~"DEG_ChromTRAP_only",
                                    TRUE~"No_DEG")

res_df$Class_selected <-case_when(res_df$Anno_threeColor == "DEG_FACS_only"~rownames(res_df),
                                  TRUE~"")

# visualization ----------------------------------------------------------------
# volcano
g_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Anno_threeColor, label = Class_selected, alpha= Anno_threeColor)) + 
  geom_point() +
  geom_text_repel(max.overlaps = Inf) +
  theme_cowplot() +
  geom_vline(xintercept = 1, linetype = 3) +
  geom_hline(yintercept = -log10(0.07), linetype = 3) +
  
  scale_color_manual(values = c("DEG_both" = "red",
                                "DEG_FACS_only" = "blue",
                                "DEG_ChromTRAP_only" = "deeppink",
                                "No_DEG" = "gray")) +
  scale_alpha_manual(values = c("DEG_both" = 0.5,
                                "DEG_FACS_only" = 0.5,
                                "DEG_ChromTRAP_only" = 0.5,
                                "No_DEG" = 0.1)) +
  labs(subtitle = paste0("FOS FACS: ", length(which(SO_FACS_HCLow@meta.data$Sample == "Cond4h_FACS")), " cells\n",
                         "HC Low rep 1: ", length(which(SO_FACS_HCLow@meta.data$Sample == "HC_OptiPrep_1")), " cells\n",
                         "HC Low rep 2: ", length(which(SO_FACS_HCLow@meta.data$Sample == "HC_OptiPrep_2")), " cells"))

ggsave(paste0(finalres_dir, "/", region, "_Volcano_FACS_vs_HCLow.png"),
       g_volcano, width = 6, height = 4)
ggsave(paste0(finalres_dir, "/", region, "_Volcano_FACS_vs_HCLow.pdf"),
       g_volcano, width = 6, height = 4)

# venn diagram
ARG_ChromTRAP <- DEseq2res_CondHigh_vs_HCLow_DEG
ARG_FACS <- rownames(res_df)[which(res_df$Class == "Upregulated")]

p <- euler(c(ChomTRAP = length(ARG_ChromTRAP) - length(dplyr::intersect(ARG_ChromTRAP, ARG_FACS)),
             "ChomTRAP&FACS" = length(dplyr::intersect(ARG_ChromTRAP, ARG_FACS)),
             FACS = length(ARG_FACS) - length(dplyr::intersect(ARG_ChromTRAP, ARG_FACS))), 
           shape = "ellipse") %>% plot(quantities = TRUE)


pdf(paste0(finalres_dir, "/", region, "_Venn_DEG_ChromTRAP_vs_FACS.pdf"))
p
dev.off()

# scatter 
DEseq2res_CondHigh_vs_HCLow$GeneName <- rownames(DEseq2res_CondHigh_vs_HCLow)
res_df$GeneName <- rownames(res_df)

df_ChromTRAP_DEG <- data.frame(GeneName = DEseq2res_CondHigh_vs_HCLow_DEG %>% unlist())
df_ChromTRAP_DEG <- left_join(df_ChromTRAP_DEG, DEseq2res_CondHigh_vs_HCLow[,c("GeneName", "log2FoldChange")])
colnames(df_ChromTRAP_DEG)[2] <- "FC_ChromTRAP"


df_ChromTRAP_DEG <- left_join(df_ChromTRAP_DEG, res_df[,c("GeneName", "log2FoldChange")])
colnames(df_ChromTRAP_DEG)[3] <- "FC_FACS"

#df_ChromTRAP_DEG <- df_ChromTRAP_DEG[-which(is.na(df_ChromTRAP_DEG$FC_FACS)),]

cor_chromTRAP_FACS <- cor(df_ChromTRAP_DEG$FC_ChromTRAP, df_ChromTRAP_DEG$FC_FACS, method = "pearson")

g_cor_chromTRAP_FACS <- ggplot(df_ChromTRAP_DEG, aes(x = FC_ChromTRAP, y = FC_FACS)) +
  geom_point() +
  geom_abline(linetype = 3) +
  xlim(c(-1, 9)) +
  ylim(c(-1, 9)) +
  theme_cowplot() +
  xlab("log2 FC ChromTRAP") +
  ylab("log2 FC FACS") +
  annotate("text", label=paste("r =", cor_chromTRAP_FACS %>% round(digits = 3)),
           x = 1, y = 7, size = 5)

ggsave(paste0(finalres_dir, "/", region, "_Scatter_Correlation_FC.pdf"), g_cor_chromTRAP_FACS, width = 4, height = 4)
ggsave(paste0(finalres_dir, "/", region, "_Scatter_Correlation_FC.png"), g_cor_chromTRAP_FACS, width = 4, height = 4)
