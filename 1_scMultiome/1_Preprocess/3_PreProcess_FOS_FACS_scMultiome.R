#!/usr/bin/env Rscript

# ==============================================================================
# scMultiome RNA/ATAC processing and chromVAR analysis
# ------------------------------------------------------------------------------
# This script processes FOS-FACS and whole-nucleus scMultiome data. It performs
# SoupX ambient RNA correction, creates a joint Seurat object, adds the ATAC
# assay, applies RNA/ATAC QC, calls peaks with MACS2, performs RNA-based
# clustering, assigns cell-type annotations, and runs chromVAR.
#
# Edit the CONFIG section before running. Sample order is fixed because barcode
# suffixes are used to assign sample labels. 
#
# Output: Seurat object with corrected RNA counts, ATAC peaks, annotations, and
# chromVAR motif deviation scores.
# ==============================================================================

library(SoupX)
library(Seurat)
library(dplyr)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(purrr)
library(harmony)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)

# Config------------------------------------------------------------------------
region <- "Hippo" # or "BLA"

data_dir_list <- c("path_to_FACS_Cond4h", "dir_to_WholeNuc_HC_rep1", "dir_to_WholeNuc_HC_rep2") # path to CellRanger result dir

annotation_file <- "metadata/CellTypeAnnotation_FOS-FACS-sorting_scMultiome.txt"

finalres_dir <- "path_to_final_res"

fragpath <- "path_to_Aggr_Cond4hFACS_HCOptiPrep2reps/atac_fragments.tsv.gz"

count <- Read10X_h5("path_to_Aggr_Cond4hFACS_HCOptiPrep2reps/filtered_feature_bc_matrix.h5")

Feature_names <- read.delim("path_to_/filtered_feature_bc_matrix/features.tsv",
                            header = FALSE,
                            stringsAsFactors = FALSE)

macs_path <- "path_to_macs2"

# helper function---------------------------------------------------------------
SoupX_func <- function(dir){
  
  # 1) read raw/filtered, get "Gene Expression" 
  raw_list  <- Read10X(file.path(dir, "raw_feature_bc_matrix"))
  filt_list <- Read10X(file.path(dir, "filtered_feature_bc_matrix"))
  
  rna_key <- grep("Gene Expression|RNA|Gene", names(raw_list), value = TRUE)[1]
  raw_rna  <- raw_list[[rna_key]]
  filt_rna <- filt_list[[rna_key]]
  
  sc <- SoupChannel(tod = raw_rna, toc = filt_rna)
  
  # 3) clustering
  seu_tmp <- CreateSeuratObject(filt_rna, project = "sample1")
  seu_tmp <- NormalizeData(seu_tmp)
  seu_tmp <- FindVariableFeatures(seu_tmp)
  seu_tmp <- ScaleData(seu_tmp, verbose = FALSE)
  seu_tmp <- RunPCA(seu_tmp, verbose = FALSE)
  seu_tmp <- FindNeighbors(seu_tmp, dims = 1:45)
  seu_tmp <- FindClusters(seu_tmp, resolution = 0.5)
  
  clust_vec <- setNames(as.character(Idents(seu_tmp)), colnames(seu_tmp))
  sc <- setClusters(sc, clusters = clust_vec)
  
  # 4) estimate contamination ratio
  sc  <- autoEstCont(sc, tfidf = TRUE)  
  
  # 5) adjust
  adj <- adjustCounts(sc)  
  
  return(adj)
  
}

# SoupX ambient RNA correction--------------------------------------------------
SoupX_adj_counts <- lapply(data_dir_list, SoupX_func)
colnames(SoupX_adj_counts[[2]]) <- sub("-1", "-2", colnames(SoupX_adj_counts[[2]]))
colnames(SoupX_adj_counts[[3]]) <- sub("-1", "-3", colnames(SoupX_adj_counts[[3]]))

# Create Seurat object----------------------------------------------------------
SO_list <- lapply(SoupX_adj_counts, CreateSeuratObject)
SO <- merge(SO_list[[1]], SO_list[2:3])

# create ATAC assay and add it to the object
annotation <- GetGRangesFromEnsDb(EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"

# merge ATAC to SO
SO[["ATAC"]] <- CreateChromatinAssay(
  counts = count$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

# QC ---------------------------------------------------------------------------      
# by RNA 
DefaultAssay(SO) <- "RNA"

mt_gene <- Feature_names$V2[which(sub("mt-", "", Feature_names$V2) != Feature_names$V2)]
SO@meta.data$Percent_mt <- PercentageFeatureSet(SO, features = mt_gene, assay = "RNA")

SO@meta.data$Sample <- case_when(sub(".+-", "", rownames(SO@meta.data)) == "1"~"Cond4h_FACS",
                                 sub(".+-", "", rownames(SO@meta.data)) == "2"~"HC_OptiPrep_1",
                                 sub(".+-", "", rownames(SO@meta.data)) == "3"~"HC_OptiPrep_2")

# by ATAC 
DefaultAssay(SO) <- "ATAC"
SO <- NucleosomeSignal(SO)
SO <- TSSEnrichment(SO)

SO <- subset(
  x = SO,
  subset = Percent_mt < 3 &
    nFeature_RNA > 300 & 
    nFeature_RNA < 8000 &
    nCount_ATAC < 100000 &
    nCount_ATAC > 500 &
    TSS.enrichment > 1 )

# Peak Calling -----------------------------------------------------------------
Peaks <- CallPeaks(SO, macs2.path = macs_path)

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
Peaks <- keepStandardChromosomes(Peaks, pruning.mode = "coarse")
Peaks <- subsetByOverlaps(x = Peaks, ranges = Signac::blacklist_mm10, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(SO),
  features = Peaks,
  cells = colnames(SO)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
annotation <- GetGRangesFromEnsDb(EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"

SO[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = Fragments(SO),
  annotation = annotation
)

# process RNA data -------------------------------------------------------------
DefaultAssay(SO) <- "RNA"

if(REGION == "Hippo"){
  SO@assays$RNA <- split(x = SO@assays$RNA, f = SO@meta.data$Sample)
}

# log-normalization
SO <- SCTransform(SO)

DefaultAssay(SO) <- "SCT"
SO <- RunPCA(SO, verbose = FALSE)

# perform UMAP -----------------------------------------------------------------
if(region == "BLA"){
  
  
  SO <- FindNeighbors(SO, dims = 1:43, reduction = "pca")
  SO <- FindClusters(SO, resolution = 0.3, graph.name = "SCT_snn")
  SO <- RunUMAP(object = SO, dims = 1:43, reduction.name = "umap_SCT", reduction = "pca")
  
  Dimplot <- DimPlot(SO, 
                     reduction = 'umap_SCT', 
                     label = TRUE, 
                     repel = TRUE, 
                     label.size = 5,
                     cols = DiscretePalette(n = 22, palette = "polychrome"),
                     group.by = "SCT_snn_res.0.3") 
  
  Dimplot_sample <- DimPlot(SO, 
                            reduction = 'umap_SCT', 
                            label = F, 
                            repel = F, 
                            label.size = 5,
                            split.by = "Sample",
                            group.by = "Sample") 
  
} else if (region == "Hippo"){
  
  SO <- IntegrateLayers(SO, method = HarmonyIntegration,
                        orig.reduction = "pca", new.reduction = "harmony")
  
  SO <- FindNeighbors(SO, dims = 1:45, reduction = "harmony")
  SO <- FindClusters(SO, resolution = 0.2, graph.name = "SCT_snn")
  SO <- RunUMAP(object = SO, dims = 1:45, reduction.name = "umap_harmony", reduction = "harmony")
  
  Dimplot <- DimPlot(SO,
                     reduction = 'umap_harmony',
                     label = TRUE,
                     repel = TRUE,
                     label.size = 5,
                     cols = DiscretePalette(n = 22, palette = "polychrome"),
                     group.by = "SCT_snn_res.0.2")
  
  Dimplot_sample <- DimPlot(SO, 
                            reduction = 'umap_harmony', 
                            label = F, 
                            repel = F, 
                            label.size = 5,
                            group.by = "Sample") 
  
  Dimplot_sample_split <- DimPlot(SO, 
                                  reduction = 'umap_harmony', 
                                  label = F, 
                                  repel = F, 
                                  label.size = 5,
                                  group.by = "Sample",
                                  split.by = "Sample") 
}

# Annotate clusters ------------------------------------------------------------
Annotation_text <- read.table(annotation_file, header = T)
Annotation_text <- Annotation_text[which(Annotation_text$region == region),]

cluster_to_celltype <- setNames(
  Annotation_text$CellType,
  Annotation_text$ClusterID
)

clusters <- as.character(SO@meta.data$seurat_clusters)

SO@meta.data$Annotation <- cluster_to_celltype[clusters]

# run chromVAR -----------------------------------------------------------------
pfm <- getMatrixSet(
  x = "path_to_JASPAR2024.sqlite3", # https://jaspar.elixir.no/downloads/
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
SO <- AddMotifs(
  object = SO,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

SO <- RunChromVAR(
  object = SO,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  new.assay.name = "chromvar"
)

# export data
saveRDS(SO, paste0(finalres_dir, "/", region, "_", "SO_chromVAR.rds"))
