# 1_scMultiome

Prepocessing and downstream analysis for single-cel Multiome (scRNA-seq + scATAC-seq) data.  

## Flow

### `1_Preprocess`: 
Generates annotated Seurat/Signac objects from 10x Multiome data.  
  The workflow includes I2 index trimming, Cell Ranger ARC `count`/`aggr`, RNA and ATAC QC, clustering, cell-type annotation, ATAC peak calling, and chromVAR motif scoring.

  Scripts in this folder process both whole-nuclei Multiome datasets and FOS-FACS-enriched Multiome datasets, then save annotated Seurat objects for downstream DEG, motif, peak, and regulatory analyses.

### `2_CellTypeMarkerPeaks`: 
Detects cell type-specific marker peaks from annotated scMultiome objects using `FindMarkers()` on the peak assay. These marker peak sets are used downstream for nanoCT/nanoCTF cell-type annotation.

### `3_DEG_WholeNucleiSeq`: 

Runs AP-1/FOS motif activity-based DEG analysis from whole-nuclei scMultiome RNA data. Cells are split into motif-high and motif-low groups using chromVAR score `MA1141.2`, followed by Seurat `FindMarkers()` for BLA/hippocampus and edgeR pseudobulk testing for mPFC.

The folder also includes scripts to summarize DEG results across brain regions and cell types, and to perform GO Molecular Function enrichment analysis for aggregated excitatory and inhibitory DEG sets.

### `4_DEG_to_benchmark_ChromTRAP`: 

Benchmarks ChromTRAP-derived DEGs against FOS-FACS-derived DEGs.  
The scripts define AP-1/FOS motif-high or motif-low excitatory neurons using chromVAR scores, perform pseudobulk DESeq2 analysis, and compare DEG overlap and fold-change consistency between ChromTRAP and FOS-FACS approaches.

### `5_ScenicPlus`: 

Runs SCENIC+ / pycisTopic-based regulatory analysis using scMultiome RNA and ATAC data.

This folder includes scripts to prepare SCENIC+ inputs, generate consensus peaks, create custom cisTarget databases, run pycisTopic and SCENIC+, and export eRegulon results. It also contains downstream analyses focused on AP-1-associated eRegulons, including enhancer classification, motif scanning with FIMO, and H3K27ac signal analysis over eRegulon enhancers.

Main outputs include curated AP-1 eRegulon tables, common/specific enhancer sets, motif enrichment summaries, and eRegulon-associated H3K27ac signal profiles.


### `6_Learning_associated_genes`

The script intersects excitatory-neuron upregulated DEGs with AP-1 eRegulon target genes, defines FOS/AP-1 motif-high and motif-low cells by chromVAR score, and clusters genes by their mean expression patterns across ChromTRAP groups. Genes expressed more highly in `Cond4h_High` than in `HC_High` cells are defined as learning-associated genes.

