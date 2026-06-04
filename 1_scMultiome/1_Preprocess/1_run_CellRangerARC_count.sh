#!/bin/bash

# ==============================================================================
# Run Cell Ranger ARC count
# ------------------------------------------------------------------------------
# This template script runs Cell Ranger ARC count for each sample listed in
# metadata/cellranger_arc_samples.tsv.
#
# Before running:
#   1. Create one libraries.csv file per sample using the template provided in
#      config/templates/cellranger_arc_libraries_template.csv.
#   2. Edit the fastqs column in each libraries.csv to match the local FASTQ
#      directory.
#   3. List each sample_id and libraries.csv path in
#      metadata/cellranger_arc_samples.txt.
#   4. Set REFERENCE, RESULT_DIR, and CELLRANGER_ARC_SIF below.
# ==============================================================================

SAMPLE_TABLE="metadata/cellranger_arc_samples.txt"
REFERENCE="path_to_cellranger_arc_reference"
RESULT_DIR="results/CellRanger"
CELLRANGER_ARC_SIF="$HOME/cellranger-arc_2.0.2.sif"

LOCALCORES=64
LOCALMEM=512

mkdir -p "${RESULT_DIR}"

tail -n +2 "${SAMPLE_TABLE}" | while IFS=$'\t' read -r sample_id libraries_csv
do
  echo "Running Cell Ranger ARC for ${sample_id}"
  echo "  libraries: ${libraries_csv}"

  if [[ ! -f "${libraries_csv}" ]]; then
  echo "ERROR: libraries.csv not found: ${libraries_csv}" >&2
  echo "Create it from config/templates/cellranger_arc_libraries_template.csv or update ${SAMPLE_TABLE}." >&2
  exit 1
  fi

  singularity exec "${CELLRANGER_ARC_SIF}" cellranger-arc count \
    --id="${sample_id}" \
    --reference="${REFERENCE}" \
    --libraries="${libraries_csv}" \
    --localcores "${LOCALCORES}" \
    --localmem "${LOCALMEM}"

  mkdir -p "${RESULT_DIR}/${sample_id}"
  mv "${sample_id}" "${RESULT_DIR}/${sample_id}/"
done
