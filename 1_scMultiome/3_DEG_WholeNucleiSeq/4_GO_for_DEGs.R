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
  
  df <- df[order(df$pvalue, decreasing = F),][1:11,]
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
DEG_Exc <- DEG_summary$Gene_name[grepl("Exc|DG|CA1|CA3", DEG_summary$celltype)]
DEG_Inh <- DEG_summary$Gene_name[grepl("Inh", DEG_summary$celltype)]

# perform GO analysis ----------------------------------------------------------
GO_Exc <- DEG_Exc %>% GO_func()
GO_Exc_plt <- GO_Exc %>% GO_plt_func()

GO_Inh <- DEG_Inh %>% GO_func()
GO_Inh_plt <- GO_Inh %>% GO_plt_func()

# Export -----------------------------------------------------------------------
ggsave(paste0(RESULT_DIR_FINAL, "/GO_bubble_ExcDEGs.png"), GO_Exc_plt, width = 6, height = 3.5)
ggsave(paste0(RESULT_DIR_FINAL, "/GO_bubble_InhDEGs.png"), GO_Inh_plt, width = 6, height = 3.5)

ggsave(paste0(RESULT_DIR_FINAL, "/GO_bubble_ExcDEGs.pdf"), GO_Exc_plt, width = 6, height = 3.5)
ggsave(paste0(RESULT_DIR_FINAL, "/GO_bubble_InhDEGs.pdf"), GO_Inh_plt, width = 6, height = 3.5)
