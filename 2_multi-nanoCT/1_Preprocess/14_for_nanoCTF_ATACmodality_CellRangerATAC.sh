#!/bin/bash

##############################################################################################################
# Aggregate replicates of Run CellRanger ATAC result from nanoscope.
# Before running, prepare cofig file by refering config/templates/cellranger_aggr_ATAC_libraries_template.csv.
##############################################################################################################

cd /path/to/your/working/dir

# config ---------------------------------------------------------------------
sampleid='samplename'
config='/path/to/your/config.csv'
reference="/path/to/your/refdata-cellranger-arc-mm10-2020-A-2.0.0"

# run cellranger atac aggr
cellranger-atac aggr --id=${sampleid} \
                     --csv=${config} \
                     --normalize=depth \
                     --reference=${reference}