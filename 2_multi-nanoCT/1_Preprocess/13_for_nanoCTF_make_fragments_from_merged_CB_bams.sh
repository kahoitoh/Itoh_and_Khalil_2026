#!/bin/bash

# ==============================================================================
# Generate fragment files from merged CB-tagged TF/IgG BAM files
# ==============================================================================
#
# Aim:
# Generate fragment files for TF and 1st antibody negative control libraries using
# CB-tagged BAM files derived from ATAC-valid-cell filtered reads.
#
# The CB-tagged BAM files are first merged across biological/technical replicates.
# Then, sinto is used to convert the merged BAM files into single-cell fragment
# files using the CB tags. The resulting fragment files are sorted, compressed,
# and indexed for downstream single-cell or peak-level analyses.
#
# Process:
# 1. Merge replicate CB-tagged BAM files for each TF sample.
# 2. Merge replicate CB-tagged BAM files for each matched 1st antibody negative
#    control sample.
# 3. Index the merged BAM files.
# 4. Generate fragment files from the merged CB-tagged BAM files using sinto.
# 5. Sort, bgzip-compress, and tabix-index the fragment files.
#
# Inputs:
# - <TF assay>_1/<TF assay>.cb.bam
# - <TF assay>_2/<TF assay>.cb.bam
# - <1st ab negative assay>_1/<1st ab negative assay>.cb.bam
# - <1st ab negative assay>_2/<1st ab negative assay>.cb.bam
#
# Outputs:
# - <TF assay>/<TF assay>.cb.bam
# - <TF assay>/<TF assay>.cb.bam.bai
# - <TF assay>/Fragments_<TF assay>.sorted.tsv.gz
# - <TF assay>/Fragments_<TF assay>.sorted.tsv.gz.tbi
# - <1st ab negative assay>/<1st ab negative assay>.cb.bam
# - <1st ab negative assay>/<1st ab negative assay>.cb.bam.bai
# - <1st ab negative assay>/Fragments_<1st ab negative assay>.sorted.tsv.gz
# - <1st ab negative assay>/Fragments_<1st ab negative assay>.sorted.tsv.gz.tbi
#
# Notes:
# - Input BAM files must already contain CB tags.
# - Replicate BAM files are merged before fragment generation.
# - sinto uses the CB tag to assign fragments to cell barcodes.
# - Fragment files are sorted by chromosome and start position before indexing.
#
# ==============================================================================

# config -----------------------------------------------------------------------
assay=(BLA_TF)  #specify sample name of TF.
assay2=(BLA_1st_ab_nega)  #specify sample name of 1st ab nega.
length=${#assay[@]}

# run --------------------------------------------------------------------------
for ((i=0; i<$length; i++)); do

    cd /path/to/output/preprocess
    mkdir ${assay[$i]} ${assay2[$i]}

    # Merge bam from replicates
    samtools merge ${assay[$i]}/${assay[$i]}.cb.bam ${assay[$i]}_1/${assay[$i]}.cb.bam ${assay[$i]}_2/${assay[$i]}.cb.bam
    samtools merge ${assay2[$i]}/${assay2[$i]}.cb.bam ${assay2[$i]}_1/${assay2[$i]}.cb.bam ${assay2[$i]}_2/${assay2[$i]}.cb.bam

    samtools index ${assay[$i]}/${assay[$i]}.cb.bam
    samtools index ${assay2[$i]}/${assay2[$i]}.cb.bam

    sinto fragments -b ${assay[$i]}/${assay[$i]} -f ${assay[$i]}/Fragments_${assay[$i]}.tsv
    cat ${assay[$i]}/Fragments_${assay[$i]}.tsv | sort -k1,1 -k2,2n | bgzip > ${assay[$i]}/Fragments_${assay[$i]}.sorted.tsv.gz
    tabix -p bed ${assay[$i]}/Fragments_${assay[$i]}.sorted.tsv.gz

    sinto fragments -b ${assay2[$i]}/${assay2[$i]} -f ${assay2[$i]}/Fragments_${assay2[$i]}.tsv
    cat ${assay2[$i]}/Fragments_${assay2[$i]}.tsv | sort -k1,1 -k2,2n | bgzip > ${assay2[$i]}/Fragments_${assay2[$i]}.sorted.tsv.gz
    tabix -p bed ${assay2[$i]}/Fragments_${assay2[$i]}.sorted.tsv.gz

done