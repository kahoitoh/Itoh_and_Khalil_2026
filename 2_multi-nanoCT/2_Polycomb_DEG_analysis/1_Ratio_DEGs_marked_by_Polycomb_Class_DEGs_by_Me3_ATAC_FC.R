#!/usr/bin/env Rscript

library(Seurat)
library(Signac)
library(dplyr)
library(GenomicRanges)
library(EnsDb.Mmusculus.v79)
library(rtracklayer)
library(purrr)
library(ggplot2)
library(edgeR)

# Config -----------------------------------------------------------------------
REGION <- c("BLA", "Hippo", "mPFC")[2]

RESULT_DIR_FINAL <- paste0("D:/in_vivo_nanoCT_2502/nCT_", REGION, "_baseres_rm_HC1_ATAC_Cells_Having_All_modalities")
RESULT_DIR_EPIC2 <- "D:/in_vivo_Multiome_EngramProject/Github/2_multi-nanoCT/1_Preprocess/data"
RESULT_DIR_DEGs_MO <- "D:/in_vivo_Multiome_EngramProject/Github/1_scMultiome/3_DEG_WholeNucleiSeq/data"
PATH_GTF <- "C:/Users/kahoi/Documents/kitazawalab_project/reference_annotation/Mus_musculus.GRCm38.102.chr.gtf" # https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/
RESULT_DIR_FIG <- "D:/in_vivo_Multiome_EngramProject/Github/2_multi-nanoCT/1_Preprocess/data"

SUBCellTYPE_LIST <- if(REGION == "BLA"){
  list("Exc" = c("LA", "BA"),
       "Inh" = c("Sst", "Vip"))
}else if(REGION == "Hippo"){
  list("Exc" = c("DG", "CA1", "CA3"),
       "Inh" = c("Inh"))
}
# Helper -----------------------------------------------------------------------
Check_gene_included_func <- function(genes_of_interest, df){
  
  df$DEG_detected <- ""
  for (g in genes_of_interest){
    
    
    df$gene <- paste0(df$overlap_gene_names)
    
    list_gene <- strsplit(df$gene, "\\|")
    index_incl_g <- lapply(list_gene, 
                           function(lg){
                             g %in% lg
                           }) %>% unlist() %>% which()
    if(length(index_incl_g) != 0){
      df$DEG_detected[index_incl_g] <- rep(g, length(index_incl_g))
    }
  }
  
  df$Chr <- sub("-.+", "", df$peak_id)
  df$Start <- sub("-[0-9]+$", "", df$peak_id) %>% sub(".+-", "", .) %>% as.numeric()
  df$End <- sub(".+-", "", df$peak_id)
  
  df$peak_id_start_corrected <- paste0(df$Chr, "-", df$Start -1, "-", df$End)
  
  df <- df[-which(df$DEG_detected == ""),]
  df
}

calc_TMM_norm_TPM <- function(RawCount){
  
  group <- factor(c("high", "low"))
  
  DGE <- DGEList(counts=RawCount %>% as.matrix(), group=group)
  DGE <- calcNormFactors(DGE)
  NormCount <- cpm(DGE, normalized.lib.sizes = TRUE)
  NormCount
 
}

calc_log2FC <- function(Count, deg){
  
  df_peaks_annotated_to_DEGs_list <- lapply(1:2, 
                                            function(n){
                                              
                                              CellID_High <- df_chromvar_count_High$CellID[which((df_chromvar_count_High$CellType %in% SUBCellTYPE_LIST[[n]]) & (df_chromvar_count_High$sample == "Cond"))]
                                              CellID_Low <- df_chromvar_count_Low$CellID[which((df_chromvar_count_Low$CellType %in% SUBCellTYPE_LIST[[n]]) & (df_chromvar_count_Low$sample == "Cond"))]
                                              
                                              df_High_Low <- data.frame(high = Count[,CellID_High] %>% rowMeans(),
                                                                        low = Count[,CellID_Low] %>% rowMeans())
                                              
                                              df_High_Low_norm <- calc_TMM_norm_TPM(df_High_Low)
                                              
                                              if(deg == T){
                                                df_peaks_annotated <- df_High_Low_norm[which(rownames(df_High_Low_norm) %in% paste0("chr", WholePeaks_anno_list[[n]]$peak_id)),] %>% as.data.frame()
                                              }else if(deg == F){
                                                df_peaks_annotated <- df_High_Low_norm[-which(rownames(df_High_Low_norm) %in% paste0("chr", WholePeaks_anno_list[[n]]$peak_id)),] %>% as.data.frame()
                                              }
                                              df_peaks_annotated$log2FC <- log2( (df_peaks_annotated[,"high"] + 1) / (df_peaks_annotated[,"low"] + 1) )
                                              df_peaks_annotated
                                            })
  df_peaks_annotated_to_DEGs <- rbind(df_peaks_annotated_to_DEGs_list[[1]], df_peaks_annotated_to_DEGs_list[[2]])
  df_peaks_annotated_to_DEGs
  
}

# Load data --------------------------------------------------------------------
combined_objects <- readRDS(paste0(RESULT_DIR_FINAL, "/", REGION, "_combined_SO_ChromVAR.rds"))

Epic2_called_H3K27me3_bed <- read.table(paste0(RESULT_DIR_EPIC2, "/", REGION, "_me_HC_merged10kb.bed"))
# Select peaks over 10kb 
Epic2_called_H3K27me3_bed <- Epic2_called_H3K27me3_bed[which((Epic2_called_H3K27me3_bed$V3 - Epic2_called_H3K27me3_bed$V2) > 10000),]

GTF <- readGFF(PATH_GTF)

DEG_all <- read.table(paste0(RESULT_DIR_DEGs_MO, "/", "DEGs_all_region_list_summary.txt"), sep = "\t", header = T)
DEG_REGION <- DEG_all[grep(REGION, DEG_all$celltype),]
DEG_Up <- DEG_REGION[which(DEG_REGION$Class_use == "Upregulated"),]

if(REGION != "Hippo"){
  
  DEGlist <- list(DEG_Up$Gene_name[grep("Exc", DEG_Up$celltype)],
                  DEG_Up$Gene_name[grep("Inh", DEG_Up$celltype)])
  
}else{
  
  DEGlist <- list(DEG_Up$Gene_name[grep("DG|CA1|CA3", DEG_Up$celltype)],
                  DEG_Up$Gene_name[grep("Inh", DEG_Up$celltype)])
  
}

annotation <- GetGRangesFromEnsDb(EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"

# Annotate peaks with overlapping genes ----------------------------------------
GTF_gene <- GTF[which(GTF$type == "gene"),]

genes_gr <- makeGRangesFromDataFrame(
  GTF_gene,
  seqnames.field = "seqid",
  start.field    = "start",
  end.field      = "end",
  strand.field   = "strand",
  keep.extra.columns = TRUE
)

prom <- promoters(genes_gr, upstream = 2000, downstream = 2000)

gene_plus_prom <- punion(genes_gr, prom, fill.gap = TRUE)

mcols(gene_plus_prom)$gene_id   <- mcols(genes_gr)$gene_id
mcols(gene_plus_prom)$gene_name <- mcols(genes_gr)$gene_name


H3K27me3_gr <-  GRanges(
  seqnames = Epic2_called_H3K27me3_bed$V1 %>% sub("chr", "", .),
  ranges   = IRanges(start = Epic2_called_H3K27me3_bed$V2 + 1,
                     end   = Epic2_called_H3K27me3_bed$V3)
)
mcols(H3K27me3_gr)$peak_id <- paste0(seqnames(H3K27me3_gr), "-", start(H3K27me3_gr), "-", end(H3K27me3_gr))
ov <- findOverlaps(H3K27me3_gr, gene_plus_prom, ignore.strand = TRUE)

ov_df <- data.frame(
  peak_idx   = queryHits(ov),
  peak_id    = mcols(H3K27me3_gr)$peak_id[queryHits(ov)],
  gene_id    = mcols(genes_gr)$gene_id[subjectHits(ov)],
  gene_name  = mcols(genes_gr)$gene_name[subjectHits(ov)]
) %>%
  distinct()

overlap_summary <- ov_df %>%
  group_by(peak_idx, peak_id) %>%
  summarise(
    overlap_gene_ids   = paste(unique(gene_id), collapse = "|"),
    overlap_gene_names = paste(unique(gene_name), collapse = "|"),
    n_overlap_genes    = n_distinct(gene_id),
    .groups = "drop"
  )

write.table(overlap_summary, 
            paste0(RESULT_DIR_FINAL,
                   "/", REGION, "_H3K27me_Annotated_over10kb.txt"), 
            quote = F, sep = "\t")

# Seach for DEGs covered by Polycomb -------------------------------------------
WholePeaks_anno_list <- lapply(DEGlist, Check_gene_included_func, overlap_summary)

# Plot ratio of Polycomb-covered DEGs ------------------------------------------
number_gene <- map2(WholePeaks_anno_list, DEGlist,
                    function(pd, d){
                      
                      c(DEG_non_peak = length(d) - length(which(d %in% pd$DEG_detected)),
                        DEG_peak = length(unique(pd$DEG_detected)))
                      
                    }) %>% unlist()

ratio_gene <- c(number_gene[1:2]/sum(number_gene[1:2]),
                number_gene[3:4]/sum(number_gene[3:4]))

df_plt_number <- data.frame(Number_DEG = number_gene,
                            Ratio_DEG = ratio_gene,
                            Class = rep(c("Exc", "Inh"), each = 2),
                            Class_DEG = names(number_gene))

g_ratio_DEG <- ggplot(df_plt_number, aes(x = Class, y = Ratio_DEG, fill = Class_DEG)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_fill_manual(values = c("DEG_peak" = "orange",
                               "DEG_non_peak" = "gray")) +
  ylab("Number of Genes") +
  xlab("") +
  labs(title = "Ratio of DEGs targeted by polycome")

ggsave(paste0(RESULT_DIR_FIG, "/", REGION, "Barplt_Ratio_of_DEGs_targeted_by_polycome_H3K27me3.png"),
       g_ratio_DEG, width = 5, height = 3)
ggsave(paste0(RESULT_DIR_FIG, "/", REGION, "Barplt_Ratio_of_DEGs_targeted_by_polycome_H3K27me3.pdf"),
       g_ratio_DEG, width = 5, height = 3)


# Quantify H3K27me3 levels in the epic2-called broad peaks ---------------------
peaks_gr <- GRanges(
  seqnames = Epic2_called_H3K27me3_bed$V1,
  ranges   = IRanges(start = Epic2_called_H3K27me3_bed$V2 + 1,
                     end   = Epic2_called_H3K27me3_bed$V3)
)

# Count me3 signals ------------------------------------------------------------
broadPeak_counts_me3 <- FeatureMatrix(
  fragments = Fragments(combined_objects$K27me3_me),
  features = peaks_gr,
  cells = colnames(combined_objects$K27me3_me)
)

combined_objects$K27me3_me[["broadPeak_H3K27me3"]] <- CreateChromatinAssay(
  counts = broadPeak_counts_me3,
  fragments = Fragments(combined_objects$K27me3_me),
  annotation = annotation
)
DefaultAssay(combined_objects$K27me3_me) <- "broadPeak_H3K27me3"

# Normalize counts
combined_objects$K27me3_me <- RunTFIDF(combined_objects$K27me3_me)

# Count ATACsignals ------------------------------------------------------------
broadPeak_counts_ATAC <- FeatureMatrix(
  fragments = Fragments(combined_objects$ATAC_aa),
  features = peaks_gr,
  cells = colnames(combined_objects$ATAC_aa)
)

combined_objects$ATAC_aa[["broadPeak_H3K27me3"]] <- CreateChromatinAssay(
  counts = broadPeak_counts_ATAC,
  fragments = Fragments(combined_objects$ATAC_aa),
  annotation = annotation
)
DefaultAssay(combined_objects$ATAC_aa) <- "broadPeak_H3K27me3"

# Normalize counts
combined_objects$ATAC_aa <- RunTFIDF(combined_objects$ATAC_aa)


# Export seurat object
saveRDS(combined_objects, paste0(RESULT_DIR_FINAL, "/", REGION, "_combined_SO_with_broadPeak_me3_Counts.rds"))

# Perform ChromTRAP ------------------------------------------------------------
# improt ChromVAR score from H3K27ac modality
# make sure the order of cells are the same between me3 and ac modalities!
combined_objects$K27me3_me@meta.data$ChromVARScore_H3K27ac <- combined_objects$K27ac_ac@assays$chromvar@data["MA1141.2",]

df_chromVAR_score <- combined_objects$K27me3_me@meta.data[,c("nCount_broadPeak_H3K27me3", "ChromVARScore_H3K27ac", "sample")]
df_chromVAR_score$sample <- sub("1|2", "", df_chromVAR_score$sample)
df_chromVAR_score$CellType <- combined_objects$K27me3_me@meta.data$Annotation
df_chromVAR_score$CellID <- rownames(combined_objects$K27me3_me@meta.data)

df_chromvar_count_High <- df_chromVAR_score %>% 
  group_by(., CellType, sample) %>% 
  mutate(q80 = quantile(ChromVARScore_H3K27ac, 0.8, na.rm = TRUE)) %>%
  dplyr::filter(ChromVARScore_H3K27ac >= q80)

df_chromvar_count_High$Condition <- "High"


df_chromvar_count_Low <- df_chromVAR_score %>% 
  group_by(., CellType, sample) %>% 
  mutate(q20 = quantile(ChromVARScore_H3K27ac, 0.2, na.rm = TRUE)) %>%
  dplyr::filter(ChromVARScore_H3K27ac <= q20)

df_chromvar_count_Low$Condition <- "Low"


# Get raw count of ATAC and me3 signal in TRAPed exc and inh neurons -----------
Count_ATAC <- GetAssayData(combined_objects$ATAC_aa, assay = "broadPeak_H3K27me3", slot = "counts")
Count_Me3 <- GetAssayData(combined_objects$K27me3_me, assay = "broadPeak_H3K27me3", slot = "counts")

# Calculate mean ATAC and me signal in exc and inh neurons,
# Select peaks associated with DEGs or non-DEGs in exc and inh neurons,
# Take CPM normalized value,
# Calculate log2FC, bind exc and inh neurons data. -----------------------------
calc_log2FC(Count_ATAC, T) %>% head()
calc_log2FC(Count_Me3, T) %>% head()

DEG_log2FC_ATAC_Me3 <- cbind(calc_log2FC(Count_ATAC, T), calc_log2FC(Count_Me3, T))
colnames(DEG_log2FC_ATAC_Me3) <- paste0(c(rep("ATAC_", 3), rep("Me3_", 3)), colnames(DEG_log2FC_ATAC_Me3))

DEG_Classed <- case_when(DEG_log2FC_ATAC_Me3$ATAC_log2FC > 0 & DEG_log2FC_ATAC_Me3$Me3_log2FC > 0~"Both_Up",
                         DEG_log2FC_ATAC_Me3$ATAC_log2FC > 0 & DEG_log2FC_ATAC_Me3$Me3_log2FC < 0~"Me_down_Ac_Up",
                         DEG_log2FC_ATAC_Me3$ATAC_log2FC < 0 & DEG_log2FC_ATAC_Me3$Me3_log2FC > 0~"Me_Up_Ac_Down",
                         DEG_log2FC_ATAC_Me3$ATAC_log2FC < 0 & DEG_log2FC_ATAC_Me3$Me3_log2FC < 0~"Both_Down") %>% table() %>% as.data.frame() %>% mutate(Ratio = Freq / sum(Freq),
                                                                                                                                                          DEG_nonDEG = "DEG")



nonDEG_log2FC_ATAC_Me3 <- cbind(calc_log2FC(Count_ATAC, F), calc_log2FC(Count_Me3, F))
colnames(nonDEG_log2FC_ATAC_Me3) <- paste0(c(rep("ATAC_", 3), rep("Me3_", 3)), colnames(nonDEG_log2FC_ATAC_Me3))
nonDEG_Classed <- case_when(nonDEG_log2FC_ATAC_Me3$ATAC_log2FC > 0 & nonDEG_log2FC_ATAC_Me3$Me3_log2FC > 0~"Both_Up",
                            nonDEG_log2FC_ATAC_Me3$ATAC_log2FC > 0 & nonDEG_log2FC_ATAC_Me3$Me3_log2FC < 0~"Me_down_Ac_Up",
                            nonDEG_log2FC_ATAC_Me3$ATAC_log2FC < 0 & nonDEG_log2FC_ATAC_Me3$Me3_log2FC > 0~"Me_Up_Ac_Down",
                            nonDEG_log2FC_ATAC_Me3$ATAC_log2FC < 0 & nonDEG_log2FC_ATAC_Me3$Me3_log2FC < 0~"Both_Down") %>% table() %>% as.data.frame() %>% mutate(Ratio = Freq / sum(Freq),
                                                                                                                                                                   DEG_nonDEG = "nonDEG")
Cat_df <- rbind(DEG_Classed, nonDEG_Classed)
colnames(Cat_df)[1] <- "Class"
Cat_df$DEG_nonDEG <- factor(Cat_df$DEG_nonDEG, levels = c("nonDEG", "DEG"))

# Perform χ² test
mat <- matrix(c(
  DEG_Classed$Freq[which(DEG_Classed[,1] == "Me_Up_Ac_Down")],
  sum(DEG_Classed$Freq) - DEG_Classed$Freq[which(DEG_Classed[,1] == "Me_Up_Ac_Down")],
  nonDEG_Classed$Freq[which(nonDEG_Classed[,1] == "Me_Up_Ac_Down")],
  sum(nonDEG_Classed$Freq) - nonDEG_Classed$Freq[which(nonDEG_Classed[,1] == "Me_Up_Ac_Down")]
), nrow = 2, byrow = TRUE)

rownames(mat) <- c("DEG", "nonDEG")
colnames(mat) <- c("Me_Up_Ac_Down", "Other")

mat
res <- fisher.test(mat) 

# Visualize
g <- ggplot(Cat_df, aes(x = DEG_nonDEG, y = Ratio, fill = Class)) +
  geom_bar(width = 0.5, stat = "identity") +
  theme_minimal() +
  scale_fill_manual(values = c("Both_Down" = "royalblue",
                               "Both_Up" = "mediumpurple1",
                               "Me_down_Ac_Up" = "maroon2",
                               "Me_Up_Ac_Down" = "lightblue2")) +
  ylab("Percentage") +
  xlab("") +
  labs(fill = "Class") +
  labs(subtitle = paste("chisq.test for Me_down_Ac_Up:", res$p.value))


ggsave(paste0(RESULT_DIR_FIG, "/", REGION, "Barplt_Ratio_of_ClassedDEGs_me3_ATAC_FC.png"),
       g, width = 5, height = 4)
ggsave(paste0(RESULT_DIR_FIG, "/", REGION, "Barplt_Ratio_of_ClassedDEGs_me3_ATAC_FC.png"),
       g, width = 5, height = 4)

