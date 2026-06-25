#!/bin/bash

# ==============================================================================
# Perform peakcalling using epic2
# ------------------------------------------------------------------------------
# This script uses epic2 to peak call H3K27me3 broad peak.
# Found peaks were stitched for 10kb using bedtools sort function.
#
# Inputs:
#   - <REGION>_me_HC.bam
#
# Outputs:
#   - <REGION>_me_HC_merged10kb.bed
#
# Edit the CONFIG section before running.
# ==============================================================================

# config --------------------------------------------------------------------------------------

cd /path/to/your/working/dir
dir_nanoscope_res='/path/to/your/results_nanoscope'

# perform epic2 and stitching peaks -----------------------------------------------------------
for region in BLA Hippo mPFC
do
    epic2 \
    --treatment ${region}_me_HC.bam \
    --genome mm10 \
    --output ${region}_me_HC_fdr005.bed

    bedtools sort -i ${region}_me_HC_fdr005.bed | bedtools merge -d 10000 -i - > ${region}_me_HC_merged10kb.bed
done