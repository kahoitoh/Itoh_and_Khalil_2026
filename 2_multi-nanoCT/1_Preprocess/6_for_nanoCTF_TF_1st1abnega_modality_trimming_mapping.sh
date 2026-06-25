#!/bin/bash

# ==============================================================================
# Trim and map multiprexed multi-nanoCT fastq files
# ------------------------------------------------------------------------------
# This script uses trimmomatic and bowtie2 to trim and map fastq files.
# This is needed for FOS and 1st-ab-negative samples because CellRanger pipeline
# does not work for them due to the hyper sparsity of the data.
#
# Inputs:
#   - .fastq from nanoscope results
#
# Outputs:
#   - <assay>/Mapped_L1.sam
#   - <assay>/Mapped_L2.sam
#
# Notes:
# Edit the CONFIG section before running.
# ==============================================================================

cd /path/to/output/preprocess

# config -----------------------------------------------------------------------
dir='/path/to/nanoscope/result/'
barcode=(GTACTGAC TATAGCCT CTAAGCCT CAGGACGT) #specity barcode. refer sub table.
assay=(BLA_1st_ab_nega_1 BLA_TF_1 BLA_1st_ab_nega_2 BLA_TF_2) #specify sample name.
fastq_prefix='nCT_TF_hM3Dq6_S1' # specify the name of .fastq.gz file.
bowtie2index='/path/to/bowtie2index' # Prepare bowtie2 index using BSgenome.Mmusculus.mm10

length=${#barcode[@]}

# perform triming and mapping --------------------------------------------------
for ((i=0; i<$length; i++)); do

    mkdir ${assay[$i]}

    trimmomatic \
    PE \
    -threads 4 \
    -phred33 \
    ${dir}/${assay[$i]}_${barcode[$i]}/fastq/barcode_${barcode[$i]}/${fastq_prefix}_L001_R1_001.fastq.gz \
    ${dir}/${assay[$i]}_${barcode[$i]}/fastq/barcode_${barcode[$i]}/${fastq_prefix}_L001_R3_001.fastq.gz \
    ${assay[$i]}/paired_L1_1_trim.fastq \
    ${assay[$i]}/unpaired_L1_1_trim.fastq \
    ${assay[$i]}/paired_L1_2_trim.fastq \
    ${assay[$i]}/unpaired_L1_2_trim.fastq \
    ILLUMINACLIP:/path/to/your/trimmomatic/adapters/NexteraPE-PE.fa:2:10:10:8:TRUE \ # modify this before running
    MINLEN:30 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15

    trimmomatic \
    PE \
    -threads 4 \
    -phred33 \
    ${dir}/${assay[$i]}_${barcode[$i]}/fastq/barcode_${barcode[$i]}/${fastq_prefix}_L002_R1_001.fastq.gz \
    ${dir}/${assay[$i]}_${barcode[$i]}/fastq/barcode_${barcode[$i]}/${fastq_prefix}_L002_R3_001.fastq.gz \
    ${assay[$i]}/paired_L2_1_trim.fastq \
    ${assay[$i]}/unpaired_L2_1_trim.fastq \
    ${assay[$i]}/paired_L2_2_trim.fastq \
    ${assay[$i]}/unpaired_L2_2_trim.fastq \
    ILLUMINACLIP:/path/to/your/trimmomatic/adapters/NexteraPE-PE.fa:2:10:10:8:TRUE \
    MINLEN:30 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15

    bowtie2 -p 16 \
    -x ${bowtie2index} \
    -1 ${assay[$i]}/paired_L1_1_trim.fastq \
    -2 ${assay[$i]}/paired_L1_2_trim.fastq \
    -S ${assay[$i]}/Mapped_L1.sam

    bowtie2 -p 16 \
    -x ${bowtie2index} \
    -1 ${assay[$i]}/paired_L2_1_trim.fastq \
    -2 ${assay[$i]}/paired_L2_2_trim.fastq \
    -S ${assay[$i]}/Mapped_L2.sam

done