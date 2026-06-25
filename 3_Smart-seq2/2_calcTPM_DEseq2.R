#!/usr/bin/env Rscript

# ==============================================================================
# Bulk RNA-seq DEG analysis and clustering
# ------------------------------------------------------------------------------
# This script processes BLA bulk RNA-seq count data, calculates TPM values,
# performs DESeq2 differential expression analysis for Conditioning samples
# against HC, and clusters conditioning-induced upregulated DEGs by their
# expression patterns across HC, Conditioning, Tone, and Shock groups.
#
# Main steps:
#   1. Filter low-count genes and calculate TPM.
#   2. Retain protein-coding genes using GTF annotation.
#   3. Run DESeq2 for Cond_1h, Cond_2h, and Cond_4h versus HC.
#   4. Summarize upregulated DEGs.
#   5. Cluster DEG mean-TPM profiles and draw a heatmap.
#
# Outputs:
#   filteredTPM.txt
#   UpDEG_summary_df.txt
#   DEseq2_Cond1h_vs_HC.txt
#   DEseq2_Cond2h_vs_HC.txt
#   DEseq2_Cond4h_vs_HC.txt
#   bulk_CondDEG_Clustering.txt
#   bulk_CondDEG_Clustering_heatmap.pdf
# ==============================================================================

library(dplyr)
library(rtracklayer)
library(DESeq2)

# config -----------------------------------------------------------------------
path_to_raw_count <- "path/to/rawCount"
path_to_reference_annotation <- "path/to/reference/annotation"
path_to_final_res <- "path/to/final"

# import rawCount
# Expected input format:
# Geneid | Chr | Start | End | Strand | Length | BLA_Cond_1h_1 | BLA_Cond_1h_2 | BLA_Cond_1h_3 | ...
CountTable_raw <- read.table(
  paste0(path_to_raw_count, "/RawCountTableContainingAllSamples.txt"),
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

# Define sample/count columns.
annotation_columns <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
count_columns <- setdiff(colnames(CountTable_raw), annotation_columns)

# remove low count genes
Count_max <- apply(CountTable_raw[,count_columns], 1, max)
Count_max %>% log10() %>% hist(breaks = 100)
length(which(Count_max < 1))
CountTable_raw_filtered <- CountTable_raw[Count_max > 1, ]

# calculate TPM ###############################################################################
Cal_TPM <- function(length_vec, rawcount){
  
  count_NormByLength <- 1000 * rawcount / length_vec
  TPM <- 1000000 * count_NormByLength / sum(count_NormByLength)
  
  return(TPM)
}

TPMTable <- CountTable_raw_filtered[, 1:6]

for (i in count_columns) {
  tpm_tmp <- Cal_TPM(
    length_vec = CountTable_raw_filtered$Length,
    rawcount = CountTable_raw_filtered[, i]
  )
  TPMTable <- cbind(TPMTable, tpm_tmp)
}

colnames(TPMTable)[(length(annotation_columns) + 1):ncol(TPMTable)] <- count_columns

write.table(TPMTable, paste0(path_to_final_res, "/filteredTPM.txt"), quote = F, sep = "\t")

# Perform DEseq2 ##############################################################################
# import GTF files
GTF <- readGFF(paste0(path_to_reference_annotation, "/Mus_musculus.GRCm38.102.chr.gtf"))
GTF_selected <- GTF[,c("gene_id", "gene_name", "gene_biotype")] %>% distinct()
colnames(GTF_selected)[1] <- "Geneid"

# For this exploratory DEG screening, only protein-coding genes were retained
# before DESeq2 analysis.
Geneid_protein_coding <- GTF_selected$Geneid[which(GTF_selected$gene_biotype == "protein_coding")] %>% unique()
CountTable_DESeq <- CountTable_raw_filtered[CountTable_raw_filtered$Geneid %in% Geneid_protein_coding,]

# Prepare count matrix for DESeq2.
count_matrix <- CountTable_DESeq[, count_columns]
rownames(count_matrix) <- CountTable_DESeq$Geneid

# Set conditions to compare with HC.
Conditions <- c("BLA_Cond_1h", "BLA_Cond_2h", "BLA_Cond_4h")

DEseq_res_list <- lapply(Conditions, 
                         function(cond, count_matrix, GTF_selected){
                           
                           RawCount_BLA_selected <- count_matrix[,grep(paste0(cond, "|HC"), colnames(count_matrix))]
                           
                           # Extract experimental condition by removing replicate-number suffix.
                           # e.g. BLA_Cond_1h_1 -> BLA_Cond_1h
                           group <- data.frame(Cond = sub("_[0-9]$", "", colnames(RawCount_BLA_selected)))
                           
                           dds <- DESeqDataSetFromMatrix(countData = RawCount_BLA_selected, colData = group, design = ~ Cond)
                           dds <- DESeq(dds) 
                           
                           
                           test_res <- results(dds, contrast = c("Cond", cond, "BLA_HC"))
                           test_res_df <- test_res %>% as.data.frame()
                           test_res_df$Geneid <- rownames(test_res_df)
                           test_res_df <- left_join(test_res_df, GTF_selected[,c("Geneid", "gene_name")], by = "Geneid")
                           
                           test_res_df <- test_res_df[!is.na(test_res_df$padj), ]
                           
                           test_res_df$Class <- case_when(test_res_df$log2FoldChange > 1.5 & test_res_df$padj < 0.0001~"Upregulated",
                                                          test_res_df$log2FoldChange < -1.5 & test_res_df$padj < 0.0001~"Downregulated",
                                                          TRUE~"Unchanged")
                           
                           return(test_res_df)
                           
                         }, count_matrix, GTF_selected)
names(DEseq_res_list) <- Conditions



# summarize DEGs
DEG_name_all_up <- c()
DEG_ID_all_up <- c()
for (i in seq_along(DEseq_res_list)){
  
  DEG_name_all_up <- c(DEG_name_all_up,
                       DEseq_res_list[[i]]$gene_name[which(DEseq_res_list[[i]]$Class == "Upregulated")])
  
  DEG_ID_all_up <- c(DEG_ID_all_up,
                     DEseq_res_list[[i]]$Geneid[which(DEseq_res_list[[i]]$Class == "Upregulated")])
  
}
DEG_summary_df <- data.frame(Geneid = DEG_ID_all_up,
                             Gene_name = DEG_name_all_up) %>% distinct()

for (i in seq_along(DEseq_res_list)){
  
  DEG_summary_df[names(DEseq_res_list)[i]] <- DEG_summary_df$Geneid %in% DEseq_res_list[[i]]$Geneid[which(DEseq_res_list[[i]]$Class == "Upregulated")]
  
}

write.table(DEG_summary_df, paste0(path_to_final_res, "/UpDEG_summary_df.txt"), quote = F, sep = "\t")
write.table(DEseq_res_list$BLA_Cond_1h, paste0(path_to_final_res, "/DEseq2_Cond1h_vs_HC.txt"), quote = F, sep = "\t")
write.table(DEseq_res_list$BLA_Cond_2h, paste0(path_to_final_res, "/DEseq2_Cond2h_vs_HC.txt"), quote = F, sep = "\t")
write.table(DEseq_res_list$BLA_Cond_4h, paste0(path_to_final_res, "/DEseq2_Cond4h_vs_HC.txt"), quote = F, sep = "\t")

# Perform clustering of DEGs ###############################################################################
samples <- c("HC",
             paste0("Cond_", c("1h", "2h", "4h")),
             paste0("Tone_", c("1h", "2h", "4h")),
             paste0("Shock_", c("1h", "2h", "4h")))

# take mean among replicates
TPM_mean <- data.frame(matrix(rep(NA, length(samples) * nrow(TPMTable)), nrow = nrow(TPMTable)))

colnames(TPM_mean) <- samples
for (i in seq_along(samples)){
  TPM_mean[,i] <- TPMTable[,grep(samples[i], colnames(TPMTable))] %>% rowMeans()
}
TPM_mean$Geneid <- TPMTable$Geneid

TPM_mean_DEG <- dplyr::filter(TPM_mean, Geneid %in% DEG_summary_df$Geneid)
TPM_mean_DEG <- left_join(TPM_mean_DEG, DEG_summary_df[,c("Geneid", "Gene_name")])

# normalization
TPM_norm <- TPM_mean_DEG[,samples] %>% t() %>% scale() %>% t() %>% as.data.frame()
rownames(TPM_norm) <- paste0(TPM_mean_DEG$Gene_name, "-", 1:nrow(TPM_mean_DEG))

# perform clustering
d <- dist(TPM_norm , method = "euclidean")
hclust_res <- hclust(d, method = "ward.D2")

TPM_clustered <- TPM_norm[hclust_res$order,]
TPM_clustered$cluster_no <- cutree(hclust_res, k = 4)[hclust_res$order]
number_unique <- cutree(hclust_res, k = 4)[hclust_res$order] %>% unique()
TPM_clustered$cluster_no_ordered <- case_when(TPM_clustered$cluster_no == number_unique[2]~1,
                                              TPM_clustered$cluster_no == number_unique[1]~2,
                                              TPM_clustered$cluster_no == number_unique[3]~3,
                                              TPM_clustered$cluster_no == number_unique[4]~4) # Reorder cluster labels for visualization in the manuscript figure.
TPM_clustered$gene_name <- rownames(TPM_clustered) %>% sub("-.+", "", .)

genes_to_label <- c("Junb", "Nr4a1", "Fos", "Dusp1", "Arc",
                    "Elmo1", "Lingo1", "Ptgs2",  "Scg2", "Pcsk1",
                    "Mlip")  
idx <- which(TPM_clustered$gene_name %in% genes_to_label)

library(ComplexHeatmap)
hm <- Heatmap(TPM_clustered[,1:10] %>% as.matrix(),
        cluster_rows = F,
        cluster_columns = F,
        heatmap_legend_param = list(title = "Z score"),
        show_row_names = F,
        row_order = 1:nrow(TPM_clustered),
        row_split = TPM_clustered$cluster_no_ordered,
        column_split = factor(c("HC", rep("Cond", 3), rep("Tone", 3), rep("Shock", 3)),
                              levels = c("HC", "Cond", "Tone", "Shock")),      
        right_annotation = rowAnnotation(
          mark = anno_mark(
            at = idx,                          
            labels = TPM_clustered$gene_name[idx],       
            labels_gp = gpar(fontsize = 10),
            link_gp = gpar(lwd = 0.6)          
          )
        ))


write.table(TPM_clustered, paste0(path_to_final_res, "/bulk_CondDEG_Clustering.txt"), quote = F, sep = "\t")

pdf(paste0(path_to_final_res, "/bulk_CondDEG_Clustering_heatmap.pdf"), width = 7, height = 8)
draw(hm)
dev.off()