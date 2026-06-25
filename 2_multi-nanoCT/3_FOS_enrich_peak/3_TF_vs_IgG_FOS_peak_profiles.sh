#!/bin/bash

# ==============================================================================
# Plot TF and 1st ab negative signal profiles around FOS enriched/non-enriched peaks
# ==============================================================================
#
# Aim:
# Compare genome-wide TF and 1st antibody negative control signals around FOS
# enriched and non-enriched peak sets using ATAC-valid-cell filtered BAM files.
#
# This script is intended to visualize whether the TF library shows stronger
# enrichment than the negative control around FOS binding peak regions.
#
# Process:
# 1. Use ATAC-valid-cell filtered, duplicate-removed, MAPQ >= 30 BAM files
#    generated in the preprocessing step.
# 2. Convert TF and 1st ab negative BAM files to RPGC-normalized bigWig files
#    using deepTools bamCoverage.
# 3. Compute signal matrices centered on:
#    - FOS enriched peaks
#    - FOS non-enriched peaks
#    using computeMatrix with +/- 3 kb windows around peak centers.
# 4. Plot average signal profiles for TF and IgG/negative control using plotProfile.
#
# Inputs:
# - <TF assay>/ATACvalidCells_MappedReads_rm_dups_q30_sort.bam
# - <1st ab negative assay>/ATACvalidCells_MappedReads_rm_dups_q30_sort.bam
# - FOS enriched peak BED files
# - FOS non-enriched peak BED files
#
# Outputs:
# - <assay>/All_ATACValidCells.bigwig
# - <assay2>/All_ATACValidCells.bigwig
# - EnrichedPeaks_edgeR_<region>_FOS_mergedRep_noLibNorm_3k.tab.gz
# - non_EnrichedPeaks_edgeR_<region>_FOS_mergedRep_noLibNorm_3k.tab.gz
# - EnrichedPeaks_edgeR_<region>_FOS_mergedRep_noLibNorm_3k.pdf
# - non_EnrichedPeaks_edgeR_<region>_FOS_mergedRep_noLibNorm_3k.pdf
#
# Notes:
# - assay and assay2 arrays must be in the same order.
# - assay contains TF sample names.
# - assay2 contains the corresponding 1st antibody negative control sample names.
# - region is used to select the corresponding FOS peak BED files.
# - BAM files must already be generated before running this script.
# - This script uses RPGC normalization in bamCoverage.
#
# ==============================================================================

# config -----------------------------------------------------------------------
assay=(BLA_TF)  #specify sample name of TF.
assay2=(BLA_1st_ab_nega)  #specify sample name of 1st ab nega.
region=(BLA)
path_to_FOS_Enrich_peak=(path_to_your_FOS_Enrich_Peak)
path_to_FOS_non_Enrich_peak=(path_to_your_FOS_non_Enrich_Peak)

length=${#assay[@]}

# run
for ((i=0; i<$length; i++)); do

    cd /path/to/output/preprocess
    
    bamCoverage \
    -b ${assay[$i]}/ATACvalidCells_MappedReads_rm_dups_q30_sort.bam  \
    -p 8 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2652783500 \
    --binSize 1 \
    -o ${assay[$i]}/All_ATACValidCells.bigwig

    bamCoverage \
    -b ${assay2[$i]}/ATACvalidCells_MappedReads_rm_dups_q30_sort.bam  \
    -p 8 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2652783500 \
    --binSize 1 \
    -o ${assay2[$i]}/All_ATACValidCells.bigwig

    computeMatrix reference-point \
    -S ${assay[$i]}/All_ATACValidCells.bigwig \
       ${assay2[$i]}/All_ATACValidCells.bigwig \
    -R ${path_to_FOS_Enrich_peak}/FOS_EnrichedPeaks_edgeR_${region[$i]}_mergedRep_fixedLibSize_TMM.bed \
    --referencePoint center \
    -a 3000 -b 3000 \
    -out EnrichedPeaks_edgeR_${region[$i]}_FOS_mergedRep_noLibNorm_3k.tab.gz

    computeMatrix reference-point \
    -S ${assay[$i]}/All_ATACValidCells.bigwig \
       ${assay2[$i]}/All_ATACValidCells.bigwig \
    -R ${path_to_FOS_non_Enrich_peak}/FOS_non_EnrichedPeaks_edgeR_${region[$i]}_mergedRep_fixedLibSize_TMM.bed \
    --referencePoint center \
    -a 3000 -b 3000 \
    -out non_EnrichedPeaks_edgeR_${region[$i]}_FOS_mergedRep_noLibNorm_3k.tab.gz

    plotProfile \
    -m EnrichedPeaks_edgeR_${region[$i]}_FOS_mergedRep_noLibNorm_3k.tab.gz \
    --perGroup \
    -out EnrichedPeaks_edgeR_${region[$i]}_FOS_mergedRep_noLibNorm_3k.pdf \
    --samplesLabel "FOS" "IgG" \
    --regionsLabel "FOS binding peaks"

    plotProfile \
    -m non_EnrichedPeaks_edgeR_${region[$i]}_FOS_mergedRep_noLibNorm_3k.tab.gz \
    --perGroup \
    -out non_EnrichedPeaks_edgeR_${region[$i]}_FOS_mergedRep_noLibNorm_3k.pdf \
    --samplesLabel "FOS" "IgG" \
    --regionsLabel "FOS binding peaks"

done