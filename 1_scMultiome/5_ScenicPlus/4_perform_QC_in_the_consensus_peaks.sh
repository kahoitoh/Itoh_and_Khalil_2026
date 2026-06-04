#!/bin/bash

pycistopic tss get_tss \
    --output outs/qc/tss.bed \
    --name "mmusculus_gene_ensembl" \
    --to-chrom-source ucsc \
    --ucsc mm10

pycistopic qc \
--fragments path_to_your_CellRanger_Aggr_dir/outs/atac_fragments.tsv.gz \
--regions path_to_your/outs/consensus_peak_calling/consensus_regions.bed \
--tss outs/qc/tss.bed \ # this is output from 
--output outs/qc/10x_multiome_brain
