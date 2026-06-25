#!/bin/bash

# Create custom cistarget database

# first, make .fa
# config
REGION_BED="/path/to/outs/consensus_peak_calling/consensus_regions.bed"
GENOME_FASTA="/path/to/fasta/mm10.fa" # wget https://hgdownload.gi.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
CHROMSIZES="/path/to/chromsize/mm10.chrom.sizes" # wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes
DATABASE_PREFIX="10x_brain_1kb_bg_with_mask"
SCRIPT_DIR="/path/to/create_cisTarget_databases" # git clone https://github.com/aertslab/create_cisTarget_databases

${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh \
        ${GENOME_FASTA} \
        ${CHROMSIZES} \
        ${REGION_BED} \
        mm10.10x_brain.with_1kb_bg_padding.fa \
        1000 \
        yes

# second, create cistarget motif database
ls /path/to/aertslab_motif_collection/v10nr_clust_public/singletons > motifs.txt # refer https://scenicplus.readthedocs.io/en/latest/human_cerebellum_ctx_db.html#Creating-custom-cistarget-database

CBDIR="path_to_/v10nr_clust_public/singletons"
FASTA_FILE="mm10.10x_brain.with_1kb_bg_padding.fa"
MOTIF_LIST="motifs.txt"
CBUST="/path/to/your/SCENIC_data/cbust" # modify this before running
DATABASE_PREFIX="10x_brain_1kb_bg_with_mask"
SCRIPT_DIR="/path/to/your/SCENIC_data/create_cisTarget_databases"  # modify this before running

"${SCRIPT_DIR}/create_cistarget_motif_databases.py" \
    -f ${FASTA_FILE} \
    -M ${CBDIR} \
    -m ${MOTIF_LIST} \
    -o ${DATABASE_PREFIX} \
    --bgpadding 1000 \
    -t 20 \
    -c "$CBUST"

