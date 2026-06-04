library(edgeR)
library(limma)
library(dplyr)

# Load bulk RNA-seq count table and scMultiome-derived gene sets
RawCount <- read.table("RawCountTableContainingAllSamples.txt",
                       sep = "\t", header = TRUE, check.names = FALSE)

MO_DEG_GRN_overlap_Genes <- read.table("Clustered_GRN_CAG_intersection.txt",
                                       header = TRUE, sep = "\t")

MO_CAG <- MO_DEG_GRN_overlap_Genes$GeneName[
  MO_DEG_GRN_overlap_Genes$cluster_no_ordered == 1
]

MO_non_CAG <- MO_DEG_GRN_overlap_Genes$GeneName[
  MO_DEG_GRN_overlap_Genes$cluster_no_ordered == 3
]

# Prepare BLA count matrix
RawCount_BLA <- RawCount[, grep("Geneid|gene_name|BLA", colnames(RawCount))]
colnames(RawCount_BLA) <- sub("BLA_", "", colnames(RawCount_BLA))

sample_info <- data.frame(
  sample = colnames(RawCount_BLA)[3:ncol(RawCount_BLA)],
  condition = sub("_[0-9]$", "", colnames(RawCount_BLA)[3:ncol(RawCount_BLA)])
)

sample_info$condition <- factor(
  sample_info$condition,
  levels = c(
    "HC",
    "Cond_1h", "Cond_2h", "Cond_4h",
    "Tone_1h", "Tone_2h", "Tone_4h",
    "Shock_1h", "Shock_2h", "Shock_4h"
  )
)

rownames(sample_info) <- sample_info$sample

counts_df <- RawCount_BLA[, 3:ncol(RawCount_BLA)]

# Remove duplicated gene symbols before gene set testing
counts_df <- counts_df[!duplicated(RawCount_BLA$gene_name), ]
rownames(counts_df) <- RawCount_BLA$gene_name[!duplicated(RawCount_BLA$gene_name)]

dge <- DGEList(counts = counts_df)
keep <- filterByExpr(dge, group = sample_info$condition)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

design <- model.matrix(~ 0 + condition, data = sample_info)
colnames(design) <- levels(sample_info$condition)

v <- voom(dge, design, plot = FALSE)

cont_2h <- makeContrasts(
  Cond2h_vs_other2h = Cond_2h - (Tone_2h + Shock_2h) / 2,
  levels = design
)

idx_list <- list(
  Cluster1 = which(rownames(v$E) %in% MO_CAG),
  Cluster3 = which(rownames(v$E) %in% MO_non_CAG)
)

fry_res_2h <- fry(
  y = v,
  index = idx_list,
  design = design,
  contrast = cont_2h[, "Cond2h_vs_other2h"]
)

write.table(fry_res_2h, "Fig6j_fry_gene_set_test_results.txt",
            quote = FALSE, sep = "\t")
