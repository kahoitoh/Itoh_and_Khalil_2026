#!/bin/bash

# ==============================================================================
# Copy and rename rawdata and run nanoscope
# ------------------------------------------------------------------------------
# This script copies raw data and chance raw data.
# Nanoscope expects the antibody barcode read as R2 and the genomic read as R3.
# Therefore, the original I2 and R2 files were renamed as R2 and R3, respectively.
#
# Inputs:
#   - <${prefix_fastq}>_L001_R1_001.fastq.gz
#   - <${prefix_fastq}>_L002_R1_001.fastq.gz
#   - <${prefix_fastq}>_L001_R2_001.fastq.gz
#   - <${prefix_fastq}>_L002_R2_001.fastq.gz
#   - <${prefix_fastq}>_L001_I2_001.fastq.gz
#   - <${prefix_fastq}>_L002_I2_001.fastq.gz
#
# Outputs: demultiplexed .fastq
#   - <assay>_<barcode>/fastq/barcode_<barcode>/<${fastq_prefix}>_L001_R1_001.fastq.gz
#   - <assay>_<barcode>/fastq/barcode_<barcode>/<${fastq_prefix}>_L001_R2_001.fastq.gz
#   - <assay>_<barcode>/fastq/barcode_<barcode>/<${fastq_prefix}>_L001_R3_001.fastq.gz
#   - <assay>_<barcode>/fastq/barcode_<barcode>/<${fastq_prefix}>_L002_R1_001.fastq.gz
#   - <assay>_<barcode>/fastq/barcode_<barcode>/<${fastq_prefix}>_L002_R2_001.fastq.gz
#   - <assay>_<barcode>/fastq/barcode_<barcode>/<${fastq_prefix}>_L002_R3_001.fastq.gz
#
# Notes:
# Edit the CONFIG section before running.
# Make config file before running nanoscope. refer config/templates/config_for_nanoscope_example.yaml.
# ==============================================================================


# Nanoscope expects the antibody barcode read as R2 and the genomic read as R3.
# Therefore, the original I2 and R2 files were renamed as R2 and R3, respectively.

cd path/to/copy/rawdata

# config -------------------------------------------------------------------
dir='path/to/your/rawdata'
prefix_fastq='sc_nCT_data_S1' # modify it aligning to your fastq filename
dir_nanoscope_res='/path/to/nanoscope/result/'
dir_nanoscope_sofware='/path/to/your/nanoscope'
dir_nanoscope_config='/path/to/your/nanoscope/config' # prepare it by refering config/templates/config_for_nanoscope_example.yaml

# copy rawdata -------------------------------------------------------------
cp ${dir}/sc_nCT_data* .

# rename fastq -------------------------------------------------------------
mv ${prefix_fastq}_L001_R2_001.fastq.gz ${prefix_fastq}_L001_R3_001.fastq.gz
mv ${prefix_fastq}_L002_R2_001.fastq.gz ${prefix_fastq}_L002_R3_001.fastq.gz

mv ${prefix_fastq}_L001_I2_001.fastq.gz ${prefix_fastq}_L001_R2_001.fastq.gz
mv ${prefix_fastq}_L002_I2_001.fastq.gz ${prefix_fastq}_L002_R2_001.fastq.gz

# run Nanoscope
cd ${dir_nanoscope_res}

snakemake --snakefile ${dir_nanoscope_software}/workflow/Snakefile_preprocess.smk \
          --cores 16 \
          --jobs 100 \
          -p \
          --use-conda \
          --configfile ${dir_nanoscope_config}/config.yaml

