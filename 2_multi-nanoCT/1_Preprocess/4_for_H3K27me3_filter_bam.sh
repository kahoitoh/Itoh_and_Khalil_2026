#!/bin/bash

# ==============================================================================
# Subset nanoCT K27me3 BAM files using neuron barcode lists
# ------------------------------------------------------------------------------
# This script uses sinto filterbarcodes to subset K27me3 BAM files to annotated
# neuronal cells for each brain region and HC replicate. The neuron barcode lists
# are generated from region-specific annotated Seurat objects.
#
# Inputs:
#   - CellID/<REGION>_Barcode_Neuron_<HC1/HC2>.txt
#   - Cell Ranger possorted_bam.bam files from nanoCT K27me3 libraries
#
# Outputs:
#   - Region/replicate-specific filtered BAM files in <REGION>_<HC1/HC2>/
#
# Edit the CONFIG section before running.
# ==============================================================================

# config --------------------------------------------------------------------------------------
cd /path/to/your/working/dir
dir_nanoscope_res='/path/to/your/results_nanoscope'

# result of each brain regions is under nCT_BLA/ nCT_Hippo/ nCT_mPFC

# filter bam ----------------------------------------------------------------------------------
for region in BLA Hippo mPFC
do
    for names in HC1 HC2
    do 
        sinto filterbarcodes \
        -b ${dir_nanoscope_res}/nCT_${region}/${names}_K27me3_*/cellranger/outs/possorted_bam.bam \
        -c CellID/${region}_Barcode_Neuron_${names}.txt \
        --outdir ${region}_${names} \
        --nproc 16
    done 
done

# merge replicate, sort, and index -------------------------------------------------------------
for region in BLA Hippo mPFC
do

    samtools merge ${region}_me_HC.bam ${region}_HC1/HC.bam ${region}_HC2/HC.bam
    samtools sort -O bam -T ${region}_me_HC -o ${region}_me_HC_sort.bam ${region}_me_HC.bam
    samtools index ${region}_me_HC_sort.bam

done