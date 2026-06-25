#!/bin/bash

# ==============================================================================
# Filter bam for TF and 1st ab negative using valid CellRanger ATAC barcode from ATAC modality
# ==============================================================================
# Aim:
# Keep only mapped reads associated with high-quality ATAC cell barcodes called by
# Cell Ranger. This is useful for TF/IgG libraries where Cell Ranger cell calling is
# unreliable because of low signal.
#
# Process:
# 1. Use filtered ATAC barcodes from the matched Cell Ranger result as valid cells.
# 2. Extract cell barcodes from the nanoscope barcode FASTQ, reverse-complement them,
#    and identify read IDs matching the valid ATAC barcode list.
#    In addition to keep_ids.txt, also generate keep_ids_with_barcode.txt, which stores
#    the read ID together with the matched Cell Ranger-compatible barcode. This file is
#    not required for the filtering steps below, but is created here for downstream use
#    when generating CB-tagged BAM files.
# 3. Use these read IDs to filter assay-specific mapped SAM files.
# 4. Merge lanes, sort, remove duplicates, filter MAPQ >= 30, and index the final BAM.
#
# Inputs:
#   - Cell Ranger barcodes.tsv from the matched ATAC assay.
#   - Nanoscope barcode FASTQ.
#   - <assay>/Mapped_L1.sam and <assay>/Mapped_L2.sam.
#
# Outputs:
#   - <assay>/keep_ids.txt
#   - <assay>/ATACvalidCells_MappedReads_rm_dups_q30_sort.bam
#   - <assay>/ATACvalidCells_MappedReads_rm_dups_q30_sort.bam.bai
#
# Notes:
#   - Edit the CONFIG section before running.
#   - barcode, assay, and assay2 arrays must be in the same order.
#   - Mapped_L1.sam and Mapped_L2.sam must already exist before running this script.
# ==============================================================================

cd /path/to/output/preprocess

# config -----------------------------------------------------------------------
dir='/path/to/nanoscope/result/'
path_to_CellRanger_res='/path/to/cellranger/result/'
barcode=(GTACTGAC TATAGCCT CTAAGCCT CAGGACGT) #specity barcode. refer sub table.
assay=(BLA_1st_ab_nega_1 BLA_TF_1 BLA_1st_ab_nega_2 BLA_TF_2) #specify sample name.
assay2=(BLA_ATAC_1 BLA_ATAC_1 BLA_ATAC_2 BLA_ATAC_2) #specify sample name of ATAC correspoding to assay.

length=${#barcode[@]}

# get cell barcode -------------------------------------------------------------
reverse_complement() {
  echo "$1" | rev | tr 'ACGTNacgtn' 'TGCANtgcan'
}

for ((i=0; i<$length; i++)); do

  mkdir ${assay[$i]}

  cut -f1 -d'-' ${path_to_CellRanger_res}/${assay2[$i]}/outs/filtered_peak_bc_matrix/barcodes.tsv > barcodes_clean.txt

  awk 'function revcomp(seq,   i, rev, c, base) {
  rev = ""
  for (i = length(seq); i > 0; i--) {
    base = substr(seq, i, 1)
    c = (base == "A" ? "T" : base == "T" ? "A" : base == "C" ? "G" : base == "G" ? "C" : "N")
    rev = rev c
  }
    return rev
  }

  BEGIN {
    while ((getline < "barcodes_clean.txt") > 0) {
      barcodes[$1] = 1
    }
  }
  {
    bc = substr($3, 1, 16)
    rbc = revcomp(bc)
    if (rbc in barcodes) print $1
  }' <(zcat ${dir}/${assay[$i]}_${barcode[$i]}/fastq/barcode_${barcode[$i]}/*_L001_R2_001.fastq.gz | paste - - - -) > ${assay[$i]}/keep_ids.txt

  awk 'function revcomp(seq,   i, rev, c, base) {
  rev = ""
  for (i = length(seq); i > 0; i--) {
    base = substr(seq, i, 1)
    c = (base == "A" ? "T" : base == "T" ? "A" : base == "C" ? "G" : base == "G" ? "C" : "N")
    rev = rev c
  }
    return rev
  }

  BEGIN {
    while ((getline < "barcodes_clean.txt") > 0) {
      barcodes[$1] = 1
    }
  }
  {
    bc = substr($3, 1, 16)
    rbc = revcomp(bc) 
    if (rbc in barcodes) print $1 "\t" rbc
  }' <(zcat ${dir}/${assay[$i]}_${barcode[$i]}/fastq/barcode_${barcode[$i]}/*_L001_R2_001.fastq.gz | paste - - - -) > ${assay[$i]}/keep_ids_with_barcode.txt

done

# subtract mapped bam based on keep_ids.txt ---------------------------------------------------
for ((i=0; i<$length; i++)); do

    sed 's/^@//' ${assay[$i]}/keep_ids.txt > ${assay[$i]}/keep_ids_noat.txt
    
    samtools view -f 3 -H ${assay[$i]}/Mapped_L1.sam | grep -E "^@" > ${assay[$i]}/tmp_test_header_L1.sam # both mates mapped
    samtools view -f 3 -H ${assay[$i]}/Mapped_L2.sam | grep -E "^@" > ${assay[$i]}/tmp_test_header_L2.sam 

    samtools view -f 3 ${assay[$i]}/Mapped_L1.sam | grep -Ff ${assay[$i]}/keep_ids_noat.txt > ${assay[$i]}/tmp_test_body_L1.sam
    samtools view -f 3 ${assay[$i]}/Mapped_L2.sam | grep -Ff ${assay[$i]}/keep_ids_noat.txt > ${assay[$i]}/tmp_test_body_L2.sam

    cat ${assay[$i]}/tmp_test_header_L1.sam ${assay[$i]}/tmp_test_body_L1.sam > ${assay[$i]}/ATACvalidCells_MappedReads_L1.sam
    cat ${assay[$i]}/tmp_test_header_L2.sam ${assay[$i]}/tmp_test_body_L2.sam > ${assay[$i]}/ATACvalidCells_MappedReads_L2.sam

    samtools merge ${assay[$i]}/ATACvalidCells_MappedReads.sam ${assay[$i]}/ATACvalidCells_MappedReads_L1.sam ${assay[$i]}/ATACvalidCells_MappedReads_L2.sam

    samtools sort -O bam -T ${assay[$i]}/ATACvalidCells_MappedReads_sort -o ${assay[$i]}/ATACvalidCells_MappedReads_sort.bam ${assay[$i]}/ATACvalidCells_MappedReads.sam
    samtools index ${assay[$i]}/ATACvalidCells_MappedReads_sort.bam

    picard MarkDuplicates \
    I=${assay[$i]}/ATACvalidCells_MappedReads_sort.bam \
    O=${assay[$i]}/ATACvalidCells_MappedReads_sort_rm_dups.bam \
    M=${assay[$i]}/rm_dups_report.txt \
    REMOVE_DUPLICATES=true

    samtools index ${assay[$i]}/ATACvalidCells_MappedReads_sort_rm_dups.bam

    samtools view -h -q 30 ${assay[$i]}/ATACvalidCells_MappedReads_sort_rm_dups.bam > ${assay[$i]}/ATACvalidCells_MappedReads_rm_dups_q30.bam
    samtools sort -O bam -T Mapped_sort -o ${assay[$i]}/ATACvalidCells_MappedReads_rm_dups_q30_sort.bam ${assay[$i]}/ATACvalidCells_MappedReads_rm_dups_q30.bam
    samtools index ${assay[$i]}/ATACvalidCells_MappedReads_rm_dups_q30_sort.bam

    rm ${assay[$i]}/tmp*

done
