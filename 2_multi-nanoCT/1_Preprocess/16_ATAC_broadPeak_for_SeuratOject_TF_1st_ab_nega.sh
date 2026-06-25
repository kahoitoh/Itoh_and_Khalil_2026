#!/bin/bash

# ==============================================================================
# Make borad Peak using ATAC modality for TF and 1st ab nega seurat object
# ==============================================================================

cd /path/to/your/working/dir
path_to_ATAC_rep1_CellRanger='path/to_ATAC/rep1/CellRanger'
path_to_ATAC_rep2_CellRanger='path/to_ATAC/rep2/CellRanger'

samtools merge Merged_ATAC.bam ${path_to_ATAC_rep1_CellRanger}/outs/possorted_bam.bam ${path_to_ATAC_rep2_CellRanger}/outs/possorted_bam.bam

macs2 callpeak \
-t Merged_ATAC.bam \
-p 1e-5 \
-f BAMPE \
-g mm \
--llocal 100000 \
--keep-dup 1 \
--broad-cutoff 0.1 \
--max-gap 1000 \
--broad 
