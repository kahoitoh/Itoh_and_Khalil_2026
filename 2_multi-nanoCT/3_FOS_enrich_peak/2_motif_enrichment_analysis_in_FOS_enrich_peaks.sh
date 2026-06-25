#!/bin/bash

# Perform motif enrichment analysis for FOS-enriched peaks

cd cd /path/to/output/preprocess

for name in BLA dHPC Cortex
do
    findMotifsGenome.pl FOS_EnrichedPeaks_edgeR_${name}_mergedRep_fixedLibSize_TMM.bed mm10 FOS_EnrichedPeaks_edgeR_${name} -size given -bg FOS_non_EnrichedPeaks_edgeR_${name}_mergedRep_fixedLibSize_TMM.bed
done
