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

## Scripts

| Script | Role | Main output |
|---|---|---|
| `0_download_files.sh` | Prepares the SCENIC+ working directory, downloads mm10 blacklist/chromsizes, prepares Mallet, initializes the SCENIC+ Snakemake pipeline, and downloads genome annotation files. | `mm10.blacklist.bed`, `scplus_pipeline/`, `Mallet-202108/`, `genome_annotation.tsv`, `chromsizes.tsv` |
| `1_modify_annotation_file.py` | Fixes chromosome naming in the genome annotation file by adding the `chr` prefix for compatibility with peak and chromsize files. | Modified `genome_annotation.tsv` |
| `2_Prepare_CellID_Annotation_RNAcounts.R` | Exports cell metadata, annotation tables, RNA count matrices, gene lists, and barcode lists from annotated scMultiome Seurat objects. | `TRAPed_<REGION>_<CELLTYPE>_anno.tsv`, RNA count `.mtx`, gene `.tsv`, barcode `.tsv` files |
| `3_make_consensus_peaks.py` | Generates pseudobulk ATAC profiles, calls peaks with MACS2, removes blacklist-overlapping regions, and creates consensus peak regions for pycisTopic / SCENIC+. | `outs/consensus_peak_calling/consensus_regions.bed`, pseudobulk BED/BigWig files, MACS peak outputs |
| `4_perform_QC_in_the_consensus_peaks.sh` | Runs pycisTopic QC using the consensus peak set and ATAC fragments. | `outs/qc/tss.bed`, `outs/qc/10x_multiome_brain/` |
| `5_create_custom_cistarget_database.sh` | Creates a custom cisTarget motif database from consensus peak regions using the Aerts lab cisTarget database tools. | `mm10.10x_brain.with_1kb_bg_padding.fa`, `motifs.txt`, custom cisTarget database files with prefix `10x_brain_1kb_bg_with_mask` |
| `6_create_cisTopic_object.py` | Creates cisTopic objects from ATAC fragments, runs MALLET LDA topic modeling, selects topic models, binarizes topics, imputes accessibility, and detects AP-1-high versus AP-1-low DARs. | `outs/cistopic_obj.pkl`, `outs/model.png`, `outs/topic_heatmap_Cell_cond.png`, topic BED files, DAR BED files |
| `7_select_peaks_for_pycisTopic.py` | Collects manually selected AP-1-high-associated topic peaks and DARs after inspecting the topic heatmap. | `outs/region_sets/<CELL_TYPE>_AP1_high_selected_region_sets/` |
| `8_prepare_RNAcount_for_SCENICplus.py` | Converts RNA count matrices and metadata into an AnnData object for SCENIC+. | `rna_raw_for_SCENICplus.h5ad` |
| `9_run_ScenicPlus.sh` | Runs the SCENIC+ Snakemake workflow after editing the SCENIC+ configuration file. | SCENIC+ pipeline outputs, including `scplusmdata.h5mu`, `ctx_results.hdf5`, and `dem_results.hdf5` |
| `10_export_ScenicPlus_res.py` | Exports SCENIC+ eRegulon metadata, AUC matrices, UMAP plots, and TF-to-gene network tables. | `df_direct_e_regulon_uns.csv`, `df_extended_e_regulon_uns.csv`, `df_direct_region_based_AUC.csv`, `df_extended_region_based_AUC.csv`, `umap_eRegulonregion_cond.png`, `edges_TF2G_full.csv` |
| `11_curate_AP-1_eRegulon.R` | Curates AP-1-associated eRegulons from SCENIC+ outputs. | `AP1_eRegulon_<REGION>.csv` and AP-1-associated TF-gene / enhancer tables |
| `12_Classify_Enhancers_into_common_or_specific.R` | Classifies AP-1 eRegulon enhancers as common or cell-type-specific across neuronal subtypes. | Common/specific AP-1 enhancer tables and BED files |
| `13_explore_celltype_specific_motifs_at_basal_state.R` | Explores cell-type-specific motif activity at basal state to identify candidate TF motifs associated with AP-1 enhancer programs. | Motif score plots and candidate motif summaries |
| `14_run_FIMO.sh` | Converts cell-type-specific AP-1 enhancer BED files to FASTA and scans candidate TF motifs using FIMO. | `fasta_<REGION>/Scenicplus_AP1_specific_enhancer_<CELLTYPE>.fa`, `fimo_<REGION>/`, `fimo_<REGION>_thresh_001/` |
| `15_plot_FIMO_res.R` | Summarizes and visualizes FIMO motif scanning results across cell-type-specific AP-1 enhancer sets. | FIMO summary tables and motif occurrence/enrichment plots |
| `16_H3K27ac_level_in_eRegulon_enhancer.R` | Quantifies nanoCT H3K27ac signal over AP-1 eRegulon enhancers. The same script structure can also be used for nanoCTF FOS signal analysis by replacing the H3K27ac input with background-subtracted FOS signal, calculated as FOS signal count minus 1st-antibody-negative count. | H3K27ac or background-subtracted FOS signal tables and plots over AP-1 eRegulon enhancer regions |

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
- For nanoCTF FOS signal analysis, `16_H3K27ac_level_in_eRegulon_enhancer.R` can be reused by replacing the H3K27ac input with background-subtracted FOS signal, calculated as `FOS signal count - 1st-antibody-negative count`.