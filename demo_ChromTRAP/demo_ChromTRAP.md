demo_ChromTRAP
================
Kaho Itoh
2026-06-25

- [Demo ChromTRAP analysis](#demo-chromtrap-analysis)
- [Load libraries](#load-libraries)
- [Config](#config)
- [Load demo data](#load-demo-data)
- [Define motif-high and motif-low
  cells](#define-motif-high-and-motif-low-cells)
- [Select cells for ChromTRAP
  comparison](#select-cells-for-chromtrap-comparison)
- [Differential expression analysis](#differential-expression-analysis)

# Demo ChromTRAP analysis

This demo reproduces the key ChromTRAP-based differential expression
analysis using a downsampled Seurat object.

The demo object, `SO_demo.rds`, contains cells from the BLA region
only.  
Specifically, it includes excitatory neuron subtypes from the basal
amygdala (BA) and lateral amygdala (LA). This reduced object is provided
so that the main analysis can be run quickly without requiring the full
Seurat object.

In this example, cells are classified by AP-1 / FOS motif activity using
the chromVAR deviation score for motif `MA1141.2`. Cells in the top
quantile are defined as motif-high cells, and cells in the bottom
quantile are defined as motif-low cells within each sample group.

# Load libraries

``` r
library(Seurat)
library(dplyr)
```

# Config

``` r
percentile_high <- 0.9
percentile_low <- 0.1
```

# Load demo data

``` r
SO <- readRDS("SO_demo.rds")
```

The demo object contains the following number of cells per annotation
and sample:

``` r
table(SO@meta.data$Annotation, SO@meta.data$Sample)
```

    ##           
    ##            Cond4h   HC
    ##   BA         1043  733
    ##   LA_Chst9    513  821
    ##   LA_Rorb    1177  833

The AP-1 / FOS motif activity score is extracted from the chromVAR assay
using motif ID `MA1141.2`.

``` r
chromvar_mat <- GetAssayData(
  object = SO,
  assay = "chromvar",
  slot = "data"
)

chromvar_count_selected <- chromvar_mat["MA1141.2", ]
```

# Define motif-high and motif-low cells

Cells are ranked by chromVAR score within each sample group.  
Here, the top 10% of cells are defined as motif-high cells, and the
bottom 10% are defined as motif-low cells.

``` r
df_chromvar_count <- data.frame(
  ChromVAR_count = as.numeric(chromvar_count_selected),
  CellID = names(chromvar_count_selected),
  HC_Cond4h = SO@meta.data[names(chromvar_count_selected), "Sample"]
)

df_High <- df_chromvar_count %>%
  group_by(HC_Cond4h) %>%
  mutate(q = quantile(ChromVAR_count, percentile_high, na.rm = TRUE)) %>%
  dplyr::filter(ChromVAR_count >= q)

df_Low <- df_chromvar_count %>%
  group_by(HC_Cond4h) %>%
  mutate(q = quantile(ChromVAR_count, percentile_low, na.rm = TRUE)) %>%
  dplyr::filter(ChromVAR_count <= q)
```

# Select cells for ChromTRAP comparison

In this demo, we compare Cond4h motif-high cells against HC motif-low
cells.

``` r
Cond_High_CellID <- df_High$CellID[df_High$HC_Cond4h == "Cond4h"]
HC_Low_CellID <- df_Low$CellID[df_Low$HC_Cond4h == "HC"]

cells_use <- c(Cond_High_CellID, HC_Low_CellID)

SO_trapped <- subset(SO, cells = cells_use)

SO_trapped$ChromTRAP_group <- case_when(
  colnames(SO_trapped) %in% Cond_High_CellID ~ "Cond4h_High",
  colnames(SO_trapped) %in% HC_Low_CellID ~ "HC_Low"
)

table(SO_trapped$ChromTRAP_group)
```

    ## 
    ## Cond4h_High      HC_Low 
    ##         274         239

# Differential expression analysis

Differential expression is performed between Cond4h motif-high cells and
HC motif-low cells using the SCT assay.

``` r
FM_res <- FindMarkers(
  object = SO_trapped,
  ident.1 = Cond_High_CellID,
  ident.2 = HC_Low_CellID,
  test.use = "wilcox",
  min.pct = 0.01,
  assay = "SCT"
)

head(FM_res)
```

    ##               p_val avg_log2FC pct.1 pct.2    p_val_adj
    ## Sorcs3 3.027435e-62   4.158315 0.869 0.155 6.347622e-58
    ## Malat1 7.297444e-59   1.024503 1.000 1.000 1.530055e-54
    ## Kcnip4 1.060556e-45  -1.148159 1.000 1.000 2.223669e-41
    ## Bdnf   6.819321e-44   5.111587 0.631 0.038 1.429807e-39
    ## Ptprd  3.193681e-41  -1.024171 1.000 1.000 6.696191e-37
    ## Etl4   3.785930e-39   2.850017 0.803 0.293 7.937960e-35
