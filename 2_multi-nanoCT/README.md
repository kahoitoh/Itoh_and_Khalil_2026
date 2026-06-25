# 2_multi-nanoCT

Preprocessing and downstream analysis for single-cell multiplexed nanoCT / nanoCTF data.

## Flow
### `1_Preprocess`:  

Preprocessing workflow for multiplexed nanoCT / nanoCTF data.

This folder contains scripts for raw FASTQ preparation, nanoscope preprocessing, histone-modality Seurat object construction, H3K27me3 broad peak calling, TF / 1st-antibody-negative read processing, ATAC-valid-cell barcode filtering, genome-browser track generation, fragment generation, Cell Ranger ATAC aggregation, and Seurat object construction for ATAC and TF/control modalities.

See `1_Preprocess/README.md` for the detailed workflow.

### `2_Polycomb_DEG_analysis`: 

Analyzes the relationship between Polycomb-associated H3K27me3 broad peaks and learning-related DEGs.


### `3_FOS_enrich_peak`: 

Identifies FOS-enriched peaks from nanoCTF TF and 1st-antibody-negative control data, then characterizes these peaks by motif enrichment and TF/control signal profiles.

