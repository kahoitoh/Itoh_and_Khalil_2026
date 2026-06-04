#!/bin/bash

cd "path/to/your/working/dir"

# Note:
#   This script is intended to be run in a cell-type-specific working directory.
#   Therefore, outputs are saved under outs/ without adding the cell type name.

# make folders
mkdir data adata_RNA outs scplus_pipeline tmp outs/qc mallet

wget https://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz
gunzip mm10.blacklist.bed.gz

# prepare for scenicplus pipeline
scenicplus init_snakemake --out_dir scplus_pipeline

wget https://github.com/mimno/Mallet/releases/download/v202108/Mallet-202108-bin.tar.gz
tar -xf Mallet-202108-bin.tar.gz

# install chromosizes and annotation file
cd scplus_pipeline/Snakemake

wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes

mv  mm10.chrom.sizes chromsizes.tsv

scenicplus prepare_data download_genome_annotations \
  --species mmusculus \
  --genome_annotation_out_fname genome_annotation.tsv \
  --chromsizes_out_fname chromsizes.tsv

# Workaround -------------------------------------------------------------------
# Some downstream SCENIC+ / pycisTopic steps expect chromsizes.tsv to be in a
# BED-like format with Chromosome, Start, and End columns. Reformat the downloaded
# chromosome size file accordingly.
awk 'BEGIN{OFS="\t"; print "Chromosome","Start","End"} {print $1,0,$2}' chromsizes.tsv > tmp.tsv # fix bed file. 
rm chromsizes.tsv
mv tmp.tsv chromsizes.tsv
