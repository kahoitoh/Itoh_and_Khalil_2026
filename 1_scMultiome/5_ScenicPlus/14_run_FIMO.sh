#!/bin/bash

# ==============================================================================
# Scan cell-type-specific AP-1-associated enhancers for TF motif occurrences
# ------------------------------------------------------------------------------
# Purpose:
#   This script scans cell-type-specific AP-1-associated enhancer sequences for
#   occurrences of candidate transcription factor motifs using FIMO.
#
# Main steps:
#   1. Select target cell types based on the brain region.
#   2. Convert cell-type-specific AP-1-associated enhancer BED files into FASTA
#      sequences using bedtools.
#   3. Load candidate TF motif files from the region-specific motif directory.
#   4. Run FIMO for each motif file and each cell type.
#   5. Save motif occurrence results into region-, cell-type-, and motif-specific
#      output directories.
#
# Inputs:
#   - bed_SpecificEnhancer_<REGION>_<CELLTYPE>.bed
#   - motif_<REGION>/*.meme
#   - Mouse reference genome FASTA
#
# Outputs:
#   - fasta_<REGION>/Scenicplus_AP1_specific_enhancer_<CELLTYPE>.fa
#   - fimo_<REGION>/<CELLTYPE>_<MOTIF_NAME>/
#
# Notes:
#   - BED files contain AP-1-associated enhancers classified as specific to each
#     cell type.
#   - Place candidate TF motif files in motif_<REGION>/ before running.
#   - FIMO is used to count motif occurrences within the enhancer sequences.
#   - The genome FASTA should match the genome build used for peak calling.
# ==============================================================================

# config ----------------------------------------------------------------------
REGION="BLA" # or "Hippo", "mPFC"

if [ "$REGION" = "BLA" ]; then
  CellTypes=("LA_Chst9" "BA" "BLA_Sst")

elif [ "$REGION" = "Hippo" ]; then
  CellTypes=("DG" "CA1" "CA3")

elif [ "$REGION" = "mPFC" ]; then
  CellTypes=("Exc" "Inh")
fi

genome_fa="/path/to/refdata-cellranger-arc-mm10-2020-A-2.0.0/fasta/genome.fa"
motif_dir="motif/motif_${REGION}" # manually curated .meme files based on ViolinPlt_ChromVARScore_Enriched_in_${celltype}.png. The meme files were downloaded from JASPAR2026.
fasta_dir="fasta_${REGION}"
fimo_dir="fimo_${REGION}"
fimo_dir_p_thresh_001="fimo_${REGION}_thresh_001"

# prepare output directories ---------------------------------------------------
mkdir -p "$fasta_dir" "$fimo_dir" "$fimo_dir_p_thresh_001"

# convert bed to fa for FIMO ---------------------------------------------------
for ct in "${CellTypes[@]}";do

    bed_file="bed_SpecificEnhancer_${REGION}_${ct}.bed"
    fasta_file="${fasta_dir}/Scenicplus_AP1_specific_enhancer_${ct}.fa"
    
    bedtools getfasta \
    -fi "$genome_fa" \
    -bed "$bed_file" \
    -name \
    -fo "$fasta_file"

done

# run FIMO --------------------------------------------------------------------
motif_files=("${motif_dir}"/*.meme)

for motif_file in "${motif_files[@]}"; do

  motif_name="$(basename "$motif_file" .meme)"

  for ct in "${CellTypes[@]}"; do

    fasta_file="${fasta_dir}/Scenicplus_AP1_specific_enhancer_${ct}.fa"
    out_dir="${fimo_dir}/${ct}_${motif_name}"
    out_dir_thresh_001="${fimo_dir_p_thresh_001}/${ct}_${motif_name}"

    fimo \
      --oc "$out_dir" \
      "$motif_file" \
      "$fasta_file"

    fimo \
      --oc "$out_dir_thresh_001" \
      --thresh 0.001 \
      "$motif_file" \
      "$fasta_file"

  done
done