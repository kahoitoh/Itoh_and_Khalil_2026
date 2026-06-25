#!/usr/bin/env Rscript

# ==============================================================================
# Summarize FIMO motif occurrences in cell-type-specific AP-1 enhancers
# ------------------------------------------------------------------------------
# This script counts enhancers with candidate TF motif hits in each cell type and
# tests motif enrichment using Fisher's exact test.
#
# Inputs:
#   - bed_SpecificEnhancer_<REGION>_<CELLTYPE>.bed
#   - fimo_<REGION>/<CELLTYPE>_<MOTIF>/fimo.tsv
#   - fimo_<REGION>_thresh_001/<CELLTYPE>_<MOTIF>/fimo.tsv
#
# Notes:
#   - Candidate TF motifs are manually selected.
#   - Motifs marked as "thresh_001" use FIMO results generated with p < 0.001.
# ==============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(ggplot2)

# config -----------------------------------------------------------------------
REGION <- c("BLA", "Hippo", "mPFC")[1]
OUT_DIR <- "path/to/your/final/res"
FIMO_RES_DIR <- "path/to/res"
BED_DIR <- "path/to/specific/enhancer/bed"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

CELLTYPES <- if (REGION == "BLA") {
  c("LA_Chst9", "BA", "BLA_Sst")
} else if (REGION == "Hippo") {
  c("DG", "CA1", "CA3")
} else if (REGION == "mPFC") {
  c("Exc", "Inh")
}

CELLTYPE_LABELS <- c(
  "LA_Chst9" = "LA_Chst9",
  "BA" = "BA",
  "BLA_Sst" = "Sst",
  "DG" = "DG",
  "CA1" = "CA1",
  "CA3" = "CA3",
  "Exc" = "Exc",
  "Inh" = "Inh"
)

MOTIF_CONFIG <- if (REGION == "BLA") {
  
  tibble::tribble(
    ~motif,     ~threshold,
    "MA0607.2", "default",     # excitatory-associated motif, bHLH (NEUROG family)
    "MA0521.3", "thresh_001"   # inhibitory-associated motif, bHLH (ASCL family)
  )
  
} else if (REGION == "Hippo") {
  
  tibble::tribble(
    ~motif,     ~threshold,
    "MA0036.4", "default",     # CA1-associated motif, GATA
    "MA1112.3", "default",     # CA3-associated motif, NR4A
    "MA0102.5", "default"      # DG-associated motif, CEBPA
  )
  
} else if (REGION == "mPFC") {
  
  tibble::tribble(
    ~motif,     ~threshold,
    "MA0802.2", "default",     # excitatory-associated motif, RFX
    "MA0607.2", "default",     # excitatory-associated motif, bHLH (NEUROG family)
    "MA1627.2", "default",     # excitatory-associated motif, WT1
    "MA0510.3", "default",     # excitatory-associated motif, TBR
    "MA0521.3", "thresh_001"   # inhibitory-associated motif, bHLH (ASCL family)
  )
  
}

# helper -----------------------------------------------------------------------
get_fimo_path <- function(motif, threshold, celltype) {
  base_dir <- ifelse(
    threshold == "thresh_001",
    paste0(FIMO_RES_DIR, "/fimo_", REGION, "_thresh_001"),
    paste0(FIMO_RES_DIR, "/fimo_", REGION)
  )
  
  paths <- c(
    file.path(base_dir, paste0(celltype, "_", motif), "fimo.txt")
  )
}

count_fimo_hits <- function(path) {
  
  x <- read.table(path, sep = "\t", header = T)
  length(unique(x$sequence_name))
  
}

count_bed <- function(celltype) {
  bed <- paste0(BED_DIR, "/bed_SpecificEnhancer_", REGION, "_", celltype, ".bed")
  nrow(read_tsv(bed, col_names = FALSE, show_col_types = FALSE))
}

plt_bubble_func <- function(fimo_res_count_list, total_enh_number_list){
  
  df <- bind_rows(lapply(names(fimo_res_count_list), function(reg){
    v <- fimo_res_count_list[[reg]]
    tibble(region = reg, motif = names(v), a = as.integer(v))
  }))
  
  N_tbl <- tibble(
    region = names(total_enh_number_list),
    N_ct   = as.integer(unlist(total_enh_number_list))
  )
  
  df <- df %>% left_join(N_tbl, by="region")
  
  
  motif_total <- df %>%
    group_by(motif) %>%
    summarise(motif_hits_all = sum(a), .groups="drop")
  
  df <- df %>% left_join(motif_total, by="motif")
  
  N_total <- sum(N_tbl$N_ct)
  
  df2 <- df %>%
    rowwise() %>%
    mutate(
      # CT = region
      b = N_ct - a,
      # not CT = other regions
      c = motif_hits_all - a,
      N_not = N_total - N_ct,
      d = N_not - c,
      
      # OR（0割・Inf回避のために擬似カウント入れてlog2ORも作る）
      OR = (a * d) / (b * c),
      log2OR = log2(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))),
      
      # Fisher（濃縮だけ見たいなら alternative="greater"）
      p = fisher.test(matrix(c(a, b, c, d), nrow=2), alternative="greater")$p.value
    ) %>%
    ungroup() %>%
    mutate(
      q = p.adjust(p, method="BH"),
      neglog10p = -log10(p),
      neglog10q = -log10(q)
    )
  
  df2$motif <- factor(df2$motif, levels = c("MA0102.5", "MA0036.4", "MA1112.3")[3:1])
  df2$region <- factor(df2$region, levels = c("DG", "CA1", "CA3"))
  
  p1 <- ggplot(df2, aes(x=region, y=motif)) +
    geom_point(aes(size=-log10(p), color=OR),  stroke = 0) +
    scale_size_continuous(range = c(1.5, 12)) +#, limits = c(0, 2.9)
    theme_minimal(base_size = 12) +
    scale_color_continuous(oob = scales::squish) +#limits = c(0, 2.7), 
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    labs(size = "-log10(p value)", color = "Odds Ratio")
  
  return(p1)
  
}

# build count lists ------------------------------------------------------------
fimo_res_count_list <- setNames(
  lapply(CELLTYPES, function(ct) {
    
    motif_counts <- purrr::map_int(seq_len(nrow(MOTIF_CONFIG)), function(i) {
      motif_i <- MOTIF_CONFIG$motif[i]
      thresh_i <- MOTIF_CONFIG$threshold[i]
      
      fimo_path <- get_fimo_path(
        motif = motif_i,
        threshold = thresh_i,
        celltype = ct
      )
      
      count_fimo_hits(fimo_path)
    })
    
    names(motif_counts) <- MOTIF_CONFIG$motif
    motif_counts
  }),
  CELLTYPES
)


total_enh_number_list <- setNames(
  lapply(CELLTYPES, function(ct) {
    count_bed(ct)
  }),
  CELLTYPES
)

# plot -------------------------------------------------------------------------

p_bubble <- plt_bubble_func(
  fimo_res_count_list = fimo_res_count_list,
  total_enh_number_list = total_enh_number_list
)

ggsave(
  filename = file.path(OUT_DIR, paste0("FIMO_motif_enrichment_bubble_", REGION, ".pdf")),
  plot = p_bubble,
  width = 4.5,
  height = 3.5
)