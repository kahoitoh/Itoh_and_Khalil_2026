#!/bin/bash

# ==============================================================================
# Aim:
# Merge biological replicate BAM files for TF and IgG/background samples, then count
# reads in Cell Ranger-called ATAC peaks for downstream TF enrichment analysis.
#
# Process:
#   1. Merge replicate BAM files for each TF sample.
#   2. Merge replicate BAM files for the matched IgG / first-antibody-negative control.
#   3. Use bedtools multicov to count TF and IgG/background reads over matched ATAC
#      peak regions.
#   4. Save a peak-by-sample count table for downstream normalization and enrichment
#      analysis.
#
# Inputs:
#   - <TF assay>_1/ATACvalidCells_MappedReads_rm_dups_q30_sort.bam
#   - <TF assay>_2/ATACvalidCells_MappedReads_rm_dups_q30_sort.bam
#   - <IgG assay>_1/ATACvalidCells_MappedReads_rm_dups_q30_sort.bam
#   - <IgG assay>_2/ATACvalidCells_MappedReads_rm_dups_q30_sort.bam
#   - <CellRanger ATAC result>/outs/peaks.bed
#
# Outputs:
#   - <TF assay>/ATACvalidCells_MappedReads_rm_dups_q30_sort.bam
#       Merged TF BAM file.
#
#   - <IgG assay>/ATACvalidCells_MappedReads_rm_dups_q30_sort.bam
#       Merged IgG/background BAM file.
#
#   - <TF assay>/counts_ATACvalidCells_<TF assay>_in_ATACpeaks_q30_rm_dups.tsv
#       Count table containing TF and IgG/background read counts in ATAC peaks.
#
# Notes:
#   - assay  = TF sample name.
#   - assay2 = matched IgG or first-antibody-negative control sample name.
#   - assay3 = matched ATAC sample name used for Cell Ranger peak calling.
#   - The arrays assay, assay2, and assay3 must be in the same order.
#   - Replicate folders are assumed to be named <assay>_1, <assay>_2,
#     <assay2>_1, and <assay2>_2.
# ==============================================================================

# config -----------------------------------------------------------------------
assay=(BLA_TF)  #specify sample name of TF.
assay2=(BLA_1st_ab_nega)  #specify sample name of 1st ab nega.
assay3=(BLA_ATAC) # specify sample name of corresponding sample of CellRanger ATAC aggr.

length=${#assay[@]}

path_to_CellRanger_aggr_res='/path/to/cellranger/result/'

cd /path/to/output/preprocess

# run --------------------------------------------------------------------------
for ((i=0; i<$length; i++)); do

    cd /path/to/output/preprocess
    mkdir ${assay[$i]} ${assay2[$i]}

    # Merge bam from replicates
    samtools merge ${assay[$i]}/ATACvalidCells_MappedReads_rm_dups_q30_sort.bam ${assay[$i]}_1/ATACvalidCells_MappedReads_rm_dups_q30_sort.bam ${assay[$i]}_2/ATACvalidCells_MappedReads_rm_dups_q30_sort.bam
    samtools merge ${assay2[$i]}/ATACvalidCells_MappedReads_rm_dups_q30_sort.bam ${assay2[$i]}_1/ATACvalidCells_MappedReads_rm_dups_q30_sort.bam ${assay2[$i]}_2/ATACvalidCells_MappedReads_rm_dups_q30_sort.bam

    samtools index ${assay[$i]}/ATACvalidCells_MappedReads_rm_dups_q30_sort.bam 
    samtools index ${assay2[$i]}/ATACvalidCells_MappedReads_rm_dups_q30_sort.bam 

    bedtools multicov \
    -bams ${assay[$i]}/ATACvalidCells_MappedReads_rm_dups_q30_sort.bam \
    ${assay2[$i]}/ATACvalidCells_MappedReads_rm_dups_q30_sort.bam \
    -bed ${path_to_CellRanger_aggr_res}/${assay3[$i]}/outs/peaks.bed \
    > ${assay[$i]}/counts_ATACvalidCells_${assay[$i]}_IgG_in_ATACpeaks_q30_rm_dups.tsv

done
