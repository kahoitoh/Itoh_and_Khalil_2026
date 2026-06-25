# 5_ScenicPlus

SCENIC+ / pycisTopic-based regulatory analysis using scMultiome RNA and ATAC data.

This folder contains scripts to prepare SCENIC+ inputs, generate consensus peaks,
create custom cisTarget databases, run pycisTopic and SCENIC+, export eRegulon
results, and perform downstream AP-1/eRegulon-focused analyses.

## Workflow overview

The workflow is organized into two parts:

1. **SCENIC+ input preparation and execution**
   - prepare genome annotation, chromsizes, blacklist, Mallet, and SCENIC+ working directories
   - export cell annotations and RNA count matrices from annotated Seurat objects
   - generate consensus ATAC peak regions
   - perform pycisTopic QC and cisTopic topic modeling
   - create custom cisTarget motif databases
   - prepare RNA input as `h5ad`
   - run SCENIC+
   - export eRegulon metadata, AUC matrices, and TF-to-gene links

2. **AP-1/eRegulon downstream analysis**
   - curate AP-1-associated eRegulons from SCENIC+ outputs
   - classify AP-1 eRegulon enhancers as common or cell-type-specific
   - identify cell-type-enriched motifs at basal state
   - scan AP-1-associated enhancers for candidate TF motifs using FIMO
   - summarize FIMO motif occurrence and enrichment
   - quantify nanoCT H3K27ac signal over AP-1 eRegulon enhancers

## Scripts

### SCENIC+ setup and input preparation

- `0_download_files.sh`  
  Prepares the SCENIC+ working directory.  
  This script creates required folders, downloads the mm10 blacklist, installs/prepares Mallet, initializes the SCENIC+ Snakemake pipeline, downloads genome annotation files, and reformats chromosome size files for downstream SCENIC+ / pycisTopic compatibility.

- `1_modify_annotation_file.py`  
  Adjusts chromosome naming in the genome annotation file.  
  This script adds the `chr` prefix to chromosome names so that genome annotation, chromsizes, and peak files use compatible chromosome naming.

- `2_Prepare_CellID_Annotation_RNAcounts.R`  
  Exports cell metadata and RNA count matrices from annotated scMultiome Seurat objects.  
  The outputs are used as RNA and cell annotation inputs for pycisTopic and SCENIC+.

- `3_make_consensus_peaks.py`  
  Generates consensus ATAC peak regions for pycisTopic / SCENIC+.  
  This script creates pseudobulk ATAC profiles, calls peaks using MACS2, removes blacklist-overlapping regions, and writes the final consensus peak BED file.

- `4_perform_QC_in_the_consensus_peaks.sh`  
  Runs pycisTopic QC using the consensus peak set and ATAC fragment file.  
  This step calculates barcode-level QC metrics used later when creating cisTopic objects.

- `5_create_custom_cistarget_database.sh`  
  Creates a custom cisTarget motif database from the consensus peak regions.  
  This script converts consensus regions to FASTA with background padding and builds a motif ranking database using the Aerts lab cisTarget database tools.

### pycisTopic and SCENIC+

- `6_create_cisTopic_object.py`  
  Runs pycisTopic analysis for AP-1-high and AP-1-low cells.  
  This script creates cisTopic objects from ATAC fragments, fits LDA topic models with MALLET, selects a topic model, binarizes topics, imputes accessibility, identifies AP-1-high versus AP-1-low DARs, and exports topic/DAR BED files.

- `7_select_peaks_for_pycisTopic.py`  
  Collects AP-1-high-associated topic peaks and DARs.  
  After manual inspection of the topic heatmap, this script copies selected AP-1-high topic peak sets and DAR BED files into a curated output directory.

- `8_prepare_RNAcount_for_SCENICplus.py`  
  Converts RNA count matrices and metadata into an `h5ad` object for SCENIC+.  
  This script prepares normalized RNA input compatible with the SCENIC+ pipeline.

- `9_run_ScenicPlus.sh`  
  Runs the SCENIC+ Snakemake workflow.  
  The SCENIC+ configuration file should be edited before running this script.

- `10_export_ScenicPlus_res.py`  
  Exports SCENIC+ results.  
  This script writes eRegulon metadata, region-based AUC matrices, UMAP plots, and TF-to-gene network tables from the SCENIC+ output object.

### AP-1/eRegulon downstream analysis

- `11_curate_AP-1_eRegulon.R`  
  Curates AP-1-associated eRegulons from SCENIC+ outputs.  
  This script extracts AP-1-related TF-gene and TF-enhancer relationships for downstream analyses.

- `12_Classify_Enhancers_into_common_or_specific.R`  
  Classifies AP-1 eRegulon enhancers as common or cell-type-specific.  
  The resulting enhancer sets are used for motif scanning and downstream signal analyses.

- `13_explore_celltype_specific_motifs_at_basal_state.R`  
  Explores cell-type-specific motif activity at basal state.  
  This script helps identify candidate TF motifs associated with cell-type-specific AP-1 enhancer programs.

- `14_run_FIMO.sh`  
  Scans cell-type-specific AP-1-associated enhancer sequences for candidate TF motifs using FIMO.  
  This script converts enhancer BED files to FASTA and runs FIMO for manually curated motif files.

- `15_plot_FIMO_res.R`  
  Summarizes and visualizes FIMO motif scanning results.  
  This script compares motif occurrence patterns across cell-type-specific AP-1 enhancer sets.

- `16_H3K27ac_level_in_eRegulon_enhancer.R`  
  Quantifies nanoCT H3K27ac signal over AP-1 eRegulon enhancers.  
  This script evaluates H3K27ac levels at eRegulon-associated enhancer regions.

## Main outputs

- consensus peak regions for SCENIC+ / pycisTopic
- pycisTopic topic regions and AP-1-high versus AP-1-low DARs
- SCENIC+ eRegulon metadata and TF-to-gene links
- curated AP-1 eRegulon tables
- common and cell-type-specific AP-1 enhancer sets
- FIMO motif occurrence summaries
- H3K27ac signal summaries over AP-1 eRegulon enhancers

## Notes

- Several steps require manual updates to paths and configuration files before running.
- pycisTopic topic selection should be adjusted based on model evaluation and topic heatmaps.
- Candidate motif files for FIMO are manually curated and should be placed in the region-specific motif directory before running `14_run_FIMO.sh`.
- The genome build used for motif scanning and cisTarget database construction should match the genome build used for peak calling.