# 1_Preprocess

Preprocessing workflow for single-cell multiplexed nanoCT / nanoCTF data.

This folder contains scripts to process raw nanoCT/nanoCTF sequencing data, build modality-specific Seurat/Signac objects, prepare H3K27me3 broad peaks, process TF and 1st-antibody-negative libraries, generate browser tracks, and create fragment/peak-count objects for downstream analyses.

## Workflow overview

The workflow is roughly organized into four parts:

1. **Raw FASTQ and nanoscope preprocessing**
   - rename FASTQ reads into the format expected by nanoscope
   - run nanoscope preprocessing

2. **Histone-modality nanoCT preprocessing**
   - build Seurat/Signac objects for ATAC, H3K27ac, and H3K27me3 modalities
   - annotate cells using scMultiome-derived marker peaks
   - run chromVAR
   - export neuronal barcodes for H3K27me3 analyses
   - filter H3K27me3 BAM files and call broad H3K27me3 peaks

3. **nanoCTF TF / 1st-antibody-negative preprocessing**
   - trim and map TF and 1st-antibody-negative FASTQ files
   - filter mapped reads using high-quality ATAC barcodes from matched ATAC modality
   - generate smoothed bedGraph and BigWig tracks
   - count TF and control reads in Cell Ranger ATAC peaks

4. **ATAC-valid-cell fragment and Seurat object generation**
   - add CB tags to TF/control BAM files
   - generate single-cell fragment files
   - run Cell Ranger ATAC aggregation
   - build annotated ATAC and TF/control Seurat objects

## Scripts

| Script | Role | Main output |
|---|---|---|
| `1_rename_fastq_run_nanoscope.sh` | Copies raw FASTQ files, renames reads into the format expected by nanoscope, and runs the nanoscope preprocessing workflow. | Demultiplexed nanoscope FASTQ outputs |
| `2_for_histone_built_SeuratObject.R` | Builds modality-specific Seurat/Signac objects for nanoCT ATAC, H3K27ac, and H3K27me3 data, performs QC, clustering, marker-peak-based annotation, and chromVAR. | Region-specific combined nanoCT Seurat objects |
| `3_for_H3K27me3_get_barcode_HC.R` | Exports annotated neuronal cell barcodes from HC K27ac nanoCT objects for downstream H3K27me3 filtering. | `<REGION>_Barcode_Neuron_HC1.txt`, `<REGION>_Barcode_Neuron_HC2.txt` |
| `4_for_H3K27me3_filter_bam.sh` | Filters K27me3 BAM files to annotated neuronal barcodes using `sinto filterbarcodes`, then merges, sorts, and indexes HC replicate BAM files. | `<REGION>_me_HC_sort.bam` |
| `5_for_H3K27me3_run_epic2_stitch_peaks.sh` | Calls broad H3K27me3 peaks using `epic2` and stitches nearby peaks within 10 kb. | `<REGION>_me_HC_merged10kb.bed` |
| `6_for_nanoCTF_TF_1st1abnega_modality_trimming_mapping.sh` | Trims and maps sparse TF / 1st-antibody-negative libraries using Trimmomatic and Bowtie2. | `<assay>/Mapped_L1.sam`, `<assay>/Mapped_L2.sam` |
| `7_for_nanoCTF_filter_mapped_reads_by_atac_barcodes.sh.sh` | Filters mapped TF/control reads using valid Cell Ranger ATAC barcodes from the matched ATAC modality, then removes duplicates and applies MAPQ filtering. | `ATACvalidCells_MappedReads_rm_dups_q30_sort.bam` |
| `8_for_nanoCTF_make_smoothed_bedgraph.py` | Generates smoothed bedGraph tracks from filtered BAM files by counting signal around paired-end fragment ends. | `*_smoothed_window200.bedgraph` |
| `9_for_nanoCTF_convert_bedgraph_to_bigwig.sh` | Sorts and filters smoothed bedGraph files, then converts them to BigWig tracks. | `*_ATACvalidCells_q30_rm_dups_smoothed_window200.bw` |
| `10_count_tf_1st_ab_nega_reads_in_atac_peaks.sh.sh` | Merges TF/control replicate BAM files and counts reads in matched Cell Ranger ATAC peaks using `bedtools multicov`. | `counts_ATACvalidCells_<TF>_IgG_in_ATACpeaks_q30_rm_dups.tsv` |
| `11_for_nanoCTF_Add_CB_to_bam.py` | Adds CB tags to BAM records using a readID-to-cell-barcode mapping file. | CB-tagged BAM |
| `12_for_nanoCTF_run_Add_CB_to_ATACvalidCells_bam.sh` | Runs the CB-tagging script for each TF/control sample. | `<assay>.cb.bam` |
| `13_for_nanoCTF_make_fragments_from_merged_CB_bams.sh` | Merges CB-tagged replicate BAM files and generates sorted/indexed single-cell fragment files using `sinto`. | `Fragments_<assay>.sorted.tsv.gz` |
| `14_for_nanoCTF_ATACmodality_CellRangerATAC.sh` | Aggregates Cell Ranger ATAC results across replicates for the ATAC modality. | Cell Ranger ATAC `aggr` output |
| `15_for_nanoCTF_Create_SeuratObject_for_ATAC.R` | Builds an annotated scATAC Seurat/Signac object from Cell Ranger ATAC peak counts, integrates replicates, scores scMultiome-derived marker peaks, and runs chromVAR. | `<REGION>_SO_ATAC.rds` |
| `16_ATAC_broadPeak_for_SeuratOject_TF_1st_ab_nega.sh` | Calls broad ATAC peaks from merged ATAC modality BAM files for TF/control Seurat object construction. | ATAC broadPeak file |
| `17_Create_SeuratObject_for_TF_and_1st_ab_nega.R` | Builds TF and 1st-antibody-negative Seurat/Signac objects from fragment files, keeps common cells, adds the control assay, and transfers annotation/FOS motif scores from the ATAC object. | `SO_<REGION>_TF_1st_ab_nega_Annotated.rds` |

## Main outputs

- annotated nanoCT Seurat objects for histone modalities
- neuronal barcode lists for H3K27me3 BAM filtering
- stitched H3K27me3 broad peak BED files
- ATAC-valid-cell-filtered TF/control BAM files
- smoothed bedGraph and BigWig tracks
- TF/control read counts in ATAC peaks
- CB-tagged BAM files and single-cell fragment files
- annotated ATAC and TF/control Seurat objects

## Notes

- Most scripts are templates and require editing the CONFIG section before running.
- Local paths, sample names, barcode arrays, and reference genome paths should be updated for each environment.
- TF and 1st-antibody-negative libraries are filtered using valid cell barcodes from the matched ATAC modality because these sparse libraries are not suitable for direct Cell Ranger cell calling.
- Several downstream scripts assume that assay arrays and corresponding ATAC sample arrays are listed in the same order.