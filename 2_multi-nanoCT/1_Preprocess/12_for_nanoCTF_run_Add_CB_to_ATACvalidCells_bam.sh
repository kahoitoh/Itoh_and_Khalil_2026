#!/bin/bash

# ==============================================================================
# Add CB tags to ATAC-valid-cell filtered TF/IgG BAM files
# ==============================================================================
#
# Aim:
# Add cell barcode information as CB tags to BAM files that were previously
# filtered using valid ATAC cell barcodes.
#
# The filtered BAM files contain mapped reads from ATAC-valid cells, but they do
# not yet have cell barcode tags in the BAM records. This script uses the
# readID-to-barcode mapping file generated during preprocessing
# (keep_ids_with_barcode.txt) and runs 12_Add_CB_to_bam.py to attach the matched
# barcode to each read as a CB tag.
#
# Process:
# 1. Loop through TF and 1st antibody negative control samples.
# 2. For each sample, use:
#    - ATACvalidCells_MappedReads_rm_dups_q30_sort.bam as the input BAM
#    - keep_ids_with_barcode.txt as the readID-to-cell-barcode map
# 3. Run 12_Add_CB_to_bam.py to add the corresponding barcode to each read using
#    the BAM CB tag.
# 4. Write the CB-tagged BAM as <assay>.cb.bam.
#
# Inputs:
# - <assay>/ATACvalidCells_MappedReads_rm_dups_q30_sort.bam
# - <assay>/keep_ids_with_barcode.txt
# - 12_Add_CB_to_bam.py
#
# Outputs:
# - <assay>.cb.bam
#
# Notes:
# - keep_ids_with_barcode.txt must contain two columns:
#     readID<TAB>cell_barcode
# - 12_Add_CB_to_bam.py removes a leading "@" from read IDs if present, so FASTQ
#   read names and BAM query names can be matched correctly.
# - The Python script reports the number of total, tagged, and unmatched reads.
# - If downstream tools require indexed BAM files, run samtools index on the
#   output BAM files after this step.
#
# ==============================================================================

# config -----------------------------------------------------------------------
assay=(BLA_1st_ab_nega_1 BLA_TF_1 BLA_1st_ab_nega_2 BLA_TF_2) #specify sample name.
length=${#assay[@]}

cd /path/to/output/preprocess

# run --------------------------------------------------------------------------
for ((i=0; i<$length; i++)); do

    python 12_Add_CB_to_bam.py \
    -i ${assay[$i]}/ATACvalidCells_MappedReads_rm_dups_q30_sort.bam \
    -m ${assay[$i]}/keep_ids_with_barcode.txt \
    -o ${assay[$i]}/${assay[$i]}.cb.bam

done
