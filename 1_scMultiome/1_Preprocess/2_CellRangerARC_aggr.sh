#!/bin/bash

# ==============================================================================
# Run Cell Ranger ARC aggr
# ------------------------------------------------------------------------------
# This template script runs Cell Ranger ARC aggr for each aggregation listed in
# metadata/cellranger_arc_aggr_samples.txt.
#
# Before running:
# Before running:
#   1. For each row in metadata/cellranger_arc_aggr_samples.txt, create the aggr.csv
#      file listed in the aggr_csv column using the template specified in the
#      template column.
#   2. Edit local molecule_h5 paths in each aggr.csv.
#   3. Set AGGR_CSV_DIR, RESULT_DIR, and CELLRANGER_ARC_SIF below.
# ==============================================================================

AGGR_TABLE="metadata/cellranger_arc_aggr_samples.txt"
AGGR_CSV_DIR="config/cellranger_arc_aggr"
RESULT_DIR="results/CellRanger_aggr"
CELLRANGER_ARC_SIF="$HOME/cellranger-arc_2.0.2.sif"

LOCALCORES=64
LOCALMEM=512

mkdir -p "${RESULT_DIR}"

tail -n +2 "${AGGR_TABLE}" | while IFS=$'\t' read -r aggr_id aggr_csv
do
  aggr_path="${AGGR_CSV_DIR}/${aggr_csv}"

  echo "Running Cell Ranger ARC aggr for ${aggr_id}"
  echo "  aggr csv: ${aggr_path}"

  if [[ ! -f "${aggr_path}" ]]; then
    echo "ERROR: aggr.csv not found: ${aggr_path}" >&2
    echo "Create it from config/templates/cellranger_arc_aggr_template.csv and update molecule_h5 paths." >&2
    exit 1
  fi

  singularity exec "${CELLRANGER_ARC_SIF}" cellranger-arc aggr \
    --id="${aggr_id}" \
    --csv="${aggr_path}" \
    --localcores="${LOCALCORES}" \
    --localmem="${LOCALMEM}"

  mkdir -p "${RESULT_DIR}"
  mv "${aggr_id}" "${RESULT_DIR}/"
done