#!/bin/bash

# ==============================================================================
# Aim:
# Convert smoothed bedGraph files into BigWig tracks for genome-browser visualization.
#
# Process:
#   1. Enter each assay folder.
#   2. Sort the smoothed bedGraph file by chromosome and genomic position.
#   3. Keep only standard mm10 chromosomes: chr1-19, chrX, chrY, and chrM.
#   4. Convert the filtered bedGraph file to BigWig using bedGraphToBigWig.
#
# Inputs:
#   - <assay>/ATACvalidCells_MappedReads_rm_dups_q30_sort_smoothed_window200.bedgraph
#   - chrom.size.mm10
#
# Outputs:
#   - <assay>/ATACvalidCells_MappedReads_rm_dups_q30_sort_smoothed_window200_sort.bedgraph
#   - <assay>/ATACvalidCells_MappedReads_rm_dups_q30_sort_smoothed_window200_sort_filtered.bedgraph
#   - <assay>/Make_smoothed_bigwig/<assay>_ATACvalidCells_q30_rm_dups_smoothed_window200.bw
#
# Notes:
#   - Run this script from the directory containing the assay folders.
#   - chrom.size.mm10 must be available in each assay folder, or the path should be
#     changed to the correct location.
#   - The output directory Make_smoothed_bigwig must already exist, or be created
#     before running bedGraphToBigWig.
# ==============================================================================

# config -----------------------------------------------------------------------
assay=(BLA_1st_ab_nega_1 BLA_TF_1 BLA_1st_ab_nega_2 BLA_TF_2)

length=${#assay[@]}

# run --------------------------------------------------------------------------
for ((i=0; i<$length; i++)); do

    cd /path/to/output/preprocess/${assay[$i]}

    LC_ALL=C sort --parallel=8 -k1,1 -k2,2n ATACvalidCells_MappedReads_rm_dups_q30_sort_smoothed_window200.bedgraph > ATACvalidCells_MappedReads_rm_dups_q30_sort_smoothed_window200_sort.bedgraph

    grep -P "^chr([1-9]|1[0-9]|X|Y|M)\t" ATACvalidCells_MappedReads_rm_dups_q30_sort_smoothed_window200_sort.bedgraph > ATACvalidCells_MappedReads_rm_dups_q30_sort_smoothed_window200_sort_filtered.bedgraph

    bedGraphToBigWig ATACvalidCells_MappedReads_rm_dups_q30_sort_smoothed_window200_sort_filtered.bedgraph chrom.size.mm10 Make_smoothed_bigwig/${assay[$i]}_ATACvalidCells_q30_rm_dups_smoothed_window200.bw

done
