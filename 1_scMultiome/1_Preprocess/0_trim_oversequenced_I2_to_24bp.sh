#!/bin/bash

# ==============================================================================
# Trim over-sequenced I2 index reads to 24 bp
# ------------------------------------------------------------------------------
# Some I2 index FASTQ files were longer than the index length expected by
# Cell Ranger ARC because of the sequencing setup. Only the affected I2 FASTQ
# files listed in Rawdata_name_oversequenced_to_sub_24bp.txt are trimmed to the
# first 24 bp. Biological R1/R2 reads are not modified.
# ==============================================================================

INPUT_LIST="metadata/Rawdata_name_oversequenced_to_sub_24bp.txt"
FASTQ_DIR="path_to_raw_fastq"
OUTDIR="processed_fastq"
THREADS=8

mkdir -p "${OUTDIR}"

echo "Copying all FASTQ files to ${OUTDIR}"
cp "${FASTQ_DIR}"/*.fastq.gz "${OUTDIR}/"

while read -r fastq; do
  [[ -z "${fastq}" ]] && continue

  in_fastq="${FASTQ_DIR}/${fastq}"
  out_fastq="${OUTDIR}/${fastq}"

  if [[ ! -f "${in_fastq}" ]]; then
    echo "ERROR: input FASTQ not found: ${in_fastq}" >&2
    echo "If FASTQ files were renamed after download, update ${INPUT_LIST}." >&2
    exit 1
  fi

  echo "Replacing ${fastq} with 24 bp-trimmed I2 FASTQ"

  seqkit subseq \
    -j "${THREADS}" \
    -r 1:24 \
    "${in_fastq}" \
    -o "${out_fastq}"

done < "${INPUT_LIST}"

seqkit stats -T "${OUTDIR}"/*_I2_001.fastq.gz \
  > "${OUTDIR}/I2_fastq_stats_after_trimming.tsv"