#!/usr/bin/env Rscript

# ==============================================================================
# Aggregate AP-1-associated DEG counts across brain regions and cell types
# ------------------------------------------------------------------------------
# This script summarizes DEG result files from BLA, hippocampus, and mPFC analyses.
# DEG tables from Seurat FindMarkers and edgeR formats are loaded and standardized.
#
# Main steps:
#   1. Load DEG result files for each region and cell type.
#   2. Extract cell type and comparison names from file names.
#   3. Standardize DEG labels across different output formats.
#   4. Aggregate Upregulated and Downregulated gene sets across three comparisons.
#   5. Count unique upregulated and downregulated genes for each cell type.
#   6. Plot DEG counts as a signed barplot.
#
# Comparisons aggregated:
#   - CondHigh_vs_CondLow
#   - CondHigh_vs_HCLow
#   - HCHigh_vs_HCLow
#
# Notes:
#   - Upregulated genes are plotted as positive values.
#   - Downregulated genes are plotted as negative values.
#   - Genes detected in multiple comparisons are counted only once per cell type
#     and regulation direction.
#   - Hippocampal DG, CA1, and CA3 are treated as excitatory populations.
# ==============================================================================

library(dplyr)
library(ggplot2)
library(tidyverse)
library(cowplot)

# config -----------------------------------------------------------------------
REGIONS <- c("BLA", "Hippo", "mPFC")
RESULT_DIR_FINAL <- "path/to/final/res"


COMPARISONS <- c(
  "CondHigh_vs_CondLow",
  "CondHigh_vs_HCLow",
  "HCHigh_vs_HCLow"
)

CELLTYPES_TO_KEEP <- c(
  "BLA_Exc", "BLA_Inh",
  "mPFC_Exc", "mPFC_Inh",
  "Hippo_DG", "Hippo_CA1", "Hippo_CA3", "Hippo_Inh"
)

# output dir -------------------------------------------------------------------
PLOT_DIR <- file.path(RESULT_DIR_FINAL, "DEG_summary_plot")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

# helper -----------------------------------------------------------------------
parse_deg_filename <- function(file_path) {
  
  file_name <- basename(file_path) %>%
    str_remove("\\.txt$")
  
  comparison <- COMPARISONS[str_detect(file_name, COMPARISONS)]
  
  celltype <- file_name %>%
    str_remove(paste0("_?DEG_?", comparison)) %>%
    str_remove(paste0("_", comparison))
  
  tibble(
    file_path = file_path,
    file_name = file_name,
    celltype = celltype,
    comparison = comparison
  )
}

read_deg_file <- function(file_path) {
  
  meta <- parse_deg_filename(file_path)
  
  read_tsv(
    file_path,
    col_types = cols(.default = col_character())
  ) %>%
    mutate(
      file_name = meta$file_name,
      celltype = meta$celltype,
      comparison = meta$comparison
    )
}

# load DEG files ---------------------------------------------------------------
deg_files <- list.files(
  RESULT_DIR_FINAL,
  pattern = "\\.txt$",
  full.names = TRUE
)

deg_meta <- map_dfr(deg_files, parse_deg_filename) %>%
  filter(
    !is.na(comparison),
    celltype %in% CELLTYPES_TO_KEEP
  )

deg_df <- deg_meta$file_path %>%
  map_dfr(read_deg_file)

# aggregate DEG gene sets ------------------------------------------------------
deg_clean <- deg_df %>%
  mutate(
    Gene_name = as.character(Gene_name),
    
    logFC_use = case_when(
      "avg_log2FC" %in% colnames(.) & !is.na(avg_log2FC) ~ as.numeric(avg_log2FC),
      "logFC" %in% colnames(.) & !is.na(logFC) ~ as.numeric(logFC),
      TRUE ~ NA_real_
    ),
    
    padj_use = case_when(
      "p_val_adj" %in% colnames(.) & !is.na(p_val_adj) ~ as.numeric(p_val_adj),
      "FDR" %in% colnames(.) & !is.na(FDR) ~ as.numeric(FDR),
      TRUE ~ NA_real_
    ),
    
    Class_use = case_when(
      Class %in% c("Upregulated", "Downregulated") ~ Class,
      DEG %in% c("Upregulated", "Downregulated") ~ DEG,
      
      Class %in% c("Unchanged", "Unchange") ~ "Unchanged",
      DEG %in% c("Unchanged", "Unchange") ~ "Unchanged",
      
      TRUE ~ "Unchanged"
    )
  ) %>%
  filter(
    celltype %in% CELLTYPES_TO_KEEP,
    comparison %in% COMPARISONS,
    !is.na(Gene_name)
  )

deg_gene_sets <- deg_clean %>%
  filter(
    comparison %in% COMPARISONS,
    Class_use %in% c("Upregulated", "Downregulated")
  ) %>%
  distinct(celltype, Class_use, Gene_name)

write.table(deg_gene_sets, paste0(RESULT_DIR_FINAL, "/DEGs_all_region_list_summary.txt"), quote = F, sep = "\t")

# count DEG numbers ------------------------------------------------------------
deg_count_df <- deg_gene_sets %>%
  count(celltype, Class_use, name = "n_genes") %>%
  mutate(
    n_genes_signed = case_when(
      Class_use == "Upregulated" ~ n_genes,
      Class_use == "Downregulated" ~ -n_genes
    ),
    Class_use = factor(Class_use, levels = c("Upregulated", "Downregulated")),
    celltype = factor(celltype, levels = CELLTYPES_TO_KEEP),
    Exc_or_Inh = sub(".+_", "", celltype) %>% sub("CA1|CA3|DG", "Exc", .)
  )

# barplot ----------------------------------------------------------------------
p <- ggplot(deg_count_df, aes(x = celltype, y = n_genes_signed, fill = Exc_or_Inh)) +
  geom_bar(stat = "identity")  +
  scale_fill_manual(values = c("Exc" = "gray20",
                               "Inh" = "gray")) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(title = "Number of DEGs") +
  ylab("Number of DEGs") +
  geom_hline(yintercept = 0, linetype = 3)
  

ggsave(
  filename = file.path(PLOT_DIR, "DEG_count_barplot.pdf"),
  plot = p,
  width = 5,
  height = 3
)

ggsave(
  filename = file.path(PLOT_DIR, "DEG_count_barplot.png"),
  plot = p,
  width = 5,
  height = 3
)
