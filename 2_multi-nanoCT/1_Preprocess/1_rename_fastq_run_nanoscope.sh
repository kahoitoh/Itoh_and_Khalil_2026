#!/bin/bash

# Nanoscope expects the antibody barcode read as R2 and the genomic read as R3.
# Therefore, the original I2 and R2 files were renamed as R2 and R3, respectively.

dir='path_to_your_rawdata'
cd your_working_dir

cp ${dir}/sc_nCT_data* .

mv sc_nCT_data_S1_L001_R2_001.fastq.gz sc_nCT_data_S1_L001_R3_001.fastq.gz
mv sc_nCT_data_S1_L002_R2_001.fastq.gz sc_nCT_data_S1_L002_R3_001.fastq.gz

mv sc_nCT_data_S1_L001_I2_001.fastq.gz sc_nCT_data_S1_L001_R2_001.fastq.gz
mv sc_nCT_data_S1_L002_I2_001.fastq.gz sc_nCT_data_S1_L002_R2_001.fastq.gz

# run Nanoscope
cd your_dir_for_results_nanoscope

snakemake --snakefile /path_to_your_nanoscope/workflow/Snakefile_preprocess.smk \
          --cores 16 \
          --jobs 100 \
          -p \
          --use-conda \
          --configfile /path_to_your_config/config.yaml

