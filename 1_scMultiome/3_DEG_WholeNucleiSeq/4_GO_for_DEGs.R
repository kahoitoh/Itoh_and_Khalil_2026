#!/usr/bin/env Rscript

# ==============================================================================
# GO enrichment analysis for aggregated excitatory and inhibitory DEG sets
# ------------------------------------------------------------------------------
# This script performs Gene Ontology enrichment analysis for DEG sets aggregated
# across brain regions and cell types.
#
# Main steps:
#   1. Load the aggregated DEG summary table.
#   2. Separate DEG gene sets into excitatory and inhibitory groups.
#   3. Convert gene symbols to Entrez IDs.
#   4. Perform GO Molecular Function enrichment analysis using clusterProfiler.
#   5. Visualize the top enriched GO terms as bubble plots.
#   6. Save plots as PNG and PDF files.
#
# Notes:
#   - Input file: DEGs_all_region_list_summary.txt
#   - Excitatory DEG sets include Exc, DG, CA1, and CA3 populations.
#   - Inhibitory DEG sets include Inh populations.
#   - GO enrichment is performed for the Molecular Function ontology.
#   - Bubble size represents gene count per GO term.
#   - Bubble color represents -log10(p-value).
# ==============================================================================

library(ggplot2)
library(cowplot)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)

# config -----------------------------------------------------------------------
RESULT_DIR_FINAL <- "path/to/final/res"

# helper -----------------------------------------------------------------------
GO_func <- function(genename){
  
  geneid <- bitr(genename, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
  
  GO_res <- enrichGO(geneid$ENTREZID, OrgDb = "org.Mm.eg.db", pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                     ont = "MF", readable = T)
  
  GO_res_df <- GO_res@result
  
  return(GO_res_df)
  
}

GO_plt_func <- function(df){
  
  df <- df[order(df$pvalue, decreasing = F),][1:13,]
  df$Description <- factor(df$Description, levels = df$Description %>% rev())
  df$x <- "x"
  
  g <- ggplot(df, aes(x = x, y = Description, size = Count, color = -log10(pvalue))) +
    geom_point() +
    scale_color_viridis_c() +
    xlab("") +
    ylab("") +
    theme_cowplot() +
    theme(axis.text.x = element_text(color = "white"))
  
  return(g)
  
}

# load data --------------------------------------------------------------------
DEG_summary <- read.table(paste0(RESULT_DIR_FINAL, "/DEGs_all_region_list_summary.txt"))

DEG_summary_up <- DEG_summary[which(DEG_summary$Class_use == "Upregulated"),]
DEG_summary_down <- DEG_summary[which(DEG_summary$Class_use == "Downregulated"),]

DEG_Exc_up <- DEG_summary_up$Gene_name[grepl("Exc|DG|CA1|CA3", DEG_summary_up$celltype)]
DEG_Inh_up <- DEG_summary_up$Gene_name[grepl("Inh", DEG_summary_up$celltype)]
DEG_Exc_down <- DEG_summary_down$Gene_name[grepl("Exc|DG|CA1|CA3", DEG_summary_down$celltype)]
DEG_Inh_down <- DEG_summary_down$Gene_name[grepl("Inh", DEG_summary_down$celltype)]

# perform GO analysis ----------------------------------------------------------
GO_Exc_up <- DEG_Exc_up %>% GO_func()
GO_Exc_plt_up <- GO_Exc_up %>% GO_plt_func()

GO_Inh_up <- DEG_Inh_up %>% GO_func()
GO_Inh_plt_up <- GO_Inh_up %>% GO_plt_func()

GO_Exc_down <- DEG_Exc_down %>% GO_func()
GO_Exc_plt_down <- GO_Exc_down %>% GO_plt_func()

GO_Inh_down <- DEG_Inh_down %>% GO_func()
GO_Inh_plt_down <- GO_Inh_down %>% GO_plt_func()

# Export -----------------------------------------------------------------------
ggsave(paste0(RESULT_DIR_FINAL, "/GO_bubble_ExcDEGs_up.png"), GO_Exc_plt_up, width = 6, height = 3.5)
ggsave(paste0(RESULT_DIR_FINAL, "/GO_bubble_InhDEGs_up.png"), GO_Inh_plt_up, width = 6, height = 3.5)

ggsave(paste0(RESULT_DIR_FINAL, "/GO_bubble_ExcDEGs_up.pdf"), GO_Exc_plt_up, width = 6, height = 3.5)
ggsave(paste0(RESULT_DIR_FINAL, "/GO_bubble_InhDEGs_up.pdf"), GO_Inh_plt_up, width = 6, height = 3.5)

ggsave(paste0(RESULT_DIR_FINAL, "/GO_bubble_ExcDEGs_down.png"), GO_Exc_plt_down, width = 8, height = 3.5)
ggsave(paste0(RESULT_DIR_FINAL, "/GO_bubble_InhDEGs_down.png"), GO_Inh_plt_down, width = 8, height = 3.5)

ggsave(paste0(RESULT_DIR_FINAL, "/GO_bubble_ExcDEGs_down.pdf"), GO_Exc_plt_down, width = 8, height = 3.5)
ggsave(paste0(RESULT_DIR_FINAL, "/GO_bubble_InhDEGs_down.pdf"), GO_Inh_plt_down, width = 8, height = 3.5)