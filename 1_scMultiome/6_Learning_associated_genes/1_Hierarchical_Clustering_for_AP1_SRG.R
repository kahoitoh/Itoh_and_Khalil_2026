

# config -----------------------------------------------------------------------
REGION <- c("BLA", "Hippo", "mPFC")[1]

PATH_TO_DEG <- "path/to/DEGs" # from 1_scMultiome/3_DEG_WholeNucleiSeq/3_summarize_DEGs.R
PATH_TO_AP1_eRegulon <- "path/to/eRegulon" # from 1_scMultiome/5_ScenicPlus/11_curate_AP-1_eRegulon.R
PATH_TO_SO_MO <- "path/to/SO/MO/wholenuclei" # from 1_scMultiome/1_Preprocess/4_PreProcess_whole_nuclei_scMultiome_each_brain_region.R

PATH_TO_FINAL <- "path/to/eRegulon"

# import data ------------------------------------------------------------------
DEG_all <- read.table(paste0(PATH_TO_DEG, "/DEGs_all_region_list_summary.txt"), sep = "\t", header = T)

# get DEGs of excitatory neurons
DEG_selected <- if(REGION == "BLA"){
  DEG_all$Gene_name[which(DEG_all$Class_use == "Upregulated" & DEG_all$celltype == "BLA_Exc")]
}else if(REGION == "Hippo"){
  DEG_all$Gene_name[which(DEG_all$Class_use == "Upregulated" & DEG_all$celltype %in% c("Hippo_CA1", "Hippo_CA3", "Hippo_DG"))]
}else if(REGION == "mPFC"){
  DEG_all$Gene_name[which(DEG_all$Class_use == "Upregulated" & DEG_all$celltype == "mPFC_Exc")]
}

eRegulon <- read.csv(paste0(PATH_TO_AP1_eRegulon, "/AP1_eRegulon_", REGION, ".csv"))
eRegulon_gene <- eRegulon$Gene %>% unique()

SO_MO_Exc <- readRDS(paste0(PATH_TO_SO_MO, "/", REGION, "_SO_list_Annotated.rds"))[[1]] # use Exc Neu object
SO_MO_Exc <- subset(SO_MO_Exc, Annotation %in% c("BA", "LA-Chst9")) # limit cell type
SO_MO_Exc@meta.data$FOS_MotifScore <- SO_MO_Exc@assays$chromvar@data["MA1141.2",]

# take intersect of DEGs and AP-1 eRegulon -------------------------------------
SRG_AP1 <- intersect(DEG_selected, eRegulon_gene)



# chromTRAP
CellID_for_TRAP_list <- lapply(c("LA-Chst9", "BA"), 
                               function(ct){
                                 
                                 tmp_meta <- SO_MO_Exc@meta.data[,c("Annotation", "Condition", "FOS_MotifScore", "CellID")]
                                 tmp_meta <- tmp_meta[which(tmp_meta$Annotation_simple %in% ct),]
                                 tmp_meta$Sample <- sub("-1|-2", "", tmp_meta$Condition)
                                 
                                 tmp_meta <- tmp_meta %>% 
                                   group_by(Sample) %>%
                                   mutate(q90 = quantile(FOS_MotifScore, 0.9, na.rm = TRUE),
                                          q10 = quantile(FOS_MotifScore, 0.1, na.rm = TRUE))
                                 
                                 tmp_meta$Class <- case_when(tmp_meta$Sample == "HC" & tmp_meta$FOS_MotifScore >= tmp_meta$q90~"HC_High",
                                                             tmp_meta$Sample == "HC" & tmp_meta$FOS_MotifScore <= tmp_meta$q10~"HC_Low",
                                                             tmp_meta$Sample == "Cond4h" & tmp_meta$FOS_MotifScore >= tmp_meta$q90~"Cond4h_High",
                                                             tmp_meta$Sample == "Cond4h" & tmp_meta$FOS_MotifScore <= tmp_meta$q10~"Cond4h_Low",
                                                             TRUE~"Others")
                                 
                                 CellID_TRAPed <- tmp_meta[,c("Class", "CellID")][which(tmp_meta$Class != "Others"),]
                                 
                                 return(CellID_TRAPed)
                               })

CellID_for_TRAP <- rbind(CellID_for_TRAP_list[[1]], CellID_for_TRAP_list[[2]])
SO_MO_Exc_TRAP <- subset(SO_MO_Exc, CellID %in% CellID_for_TRAP$CellID)

tmp_meta <- SO_MO_Exc_TRAP@meta.data
tmp_meta <- left_join(tmp_meta, CellID_for_TRAP)
SO_MO_Exc_TRAP@meta.data$Class_TRAP <- tmp_meta$Class

# Perform hierarchical clustering of AP1-SRG among four groups defined by ChromTRAP
## get expression levels of AP-1 SRG
gene_count <- GetAssayData(SO_MO_Exc_TRAP, assay = "SCT", slot = "data")
index <- which(rownames(gene_count) %in% GRN_Gene_list)
gene_count_selected <- gene_count[index,]

cell_cond <- SO_MO_Exc_TRAP@meta.data[colnames(gene_count_selected), "Class_TRAP"]
gene_mean_by_cond <- sapply(
  split(seq_along(cell_cond), cell_cond),
  function(idx) {
    Matrix::rowMeans(gene_count_selected[, idx, drop = FALSE])
  }
) %>% as.data.frame() # calcualte mean of DEG exp in each group

gene_mean_by_cond$log2FC_Cond <- log2( (gene_mean_by_cond$Cond4h_High + 0.001) / (gene_mean_by_cond$Cond4h_Low + 0.001) )
gene_mean_by_cond$log2FC_HC <- log2( (gene_mean_by_cond$HC_High + 0.001) / (gene_mean_by_cond$HC_Low + 0.001) )

df_norm <- gene_mean_by_cond[,1:4] %>% t() %>% scale() %>% t() %>% as.data.frame()

## perform clustering
d <- dist(df_norm, method = "euclidean")
hclust_res <- hclust(d, method = "ward.D2")
df_norm <- df_norm[hclust_res$order,]

df_norm$Clustered_number <- 1:nrow(df_norm)
df_norm$GeneName <- rownames(df_norm)

df_norm$cluster_no <- cutree(hclust_res, k = 3)[hclust_res$order]
number_unique <- cutree(hclust_res, k = 3)[hclust_res$order] %>% unique()
df_norm$cluster_no_ordered <- case_when(df_norm$cluster_no == number_unique[1]~1,
                                        df_norm$cluster_no == number_unique[2]~2,
                                        df_norm$cluster_no == number_unique[3]~3)

Heatmap( df_norm[,1:4] %>% as.matrix(),
         cluster_rows = F,
         cluster_columns = F,
         heatmap_legend_param = list(title = "Z score"),
         show_row_names = F,
         row_order = 1:nrow(df_norm),
         row_split =  df_norm$cluster_no_ordered,
         use_raster = F)

write.table(df_norm,
            paste0(PATH_TO_FINAL, "/Clustered_AP1_SRG_", REGION, ".txt"),
            sep = "\t", quote = F)
