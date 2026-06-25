#!/bin/bash

cd path_to_working_dir

dir="path_to_rawdata"
gtfpath="path_to_ref"
sample_list="SampleName.txt"

mkdir -p result/trimming result/STAR result/featureCount

# trimming
while read line1 <&3; do

    trimmomatic \
    PE \
    -threads 4 \
    -phred33 \
    ${dir}/${line1}_L1_1.fq.gz \
    ${dir}/${line1}_L1_2.fq.gz \
    result/trimming/paired_${line1}_L1_1_trim.fastq \
    result/trimming/unpaired_${line1}_L1_1_trim.fastq \
    result/trimming/paired_${line1}_L1_2_trim.fastq \
    result/trimming/unpaired_${line1}_L1_2_trim.fastq \
    ILLUMINACLIP:/path/to/your/trimmomatic/adapters/NexteraPE-PE.fa:2:10:10 \ # modify this before running
    MINLEN:30 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15

done 3<"SampleName.txt"

# mapping
while read line1 <&3; do

    STAR --runThreadN 16 \
     --genomeDir /path_to_STAR_index/STARindex_Mus_musculus.GRCm38.dna.primary_assembly \
     --readFilesIn result/trimming/paired_${line1}_L1_1_trim.fastq  result/trimming/paired_${line1}_L1_2_trim.fastq \
     --genomeLoad NoSharedMemory \
     --outFilterMultimapNmax 1 \
     --outFileNamePrefix result/STAR/${line1}_L1

done 3<"SampleName.txt"

# count
gtfpath='path_to_ref'
while read line1 <&3; do

    samtools sort -O sam -T result/STAR/${line1}_sort -o result/STAR/${line1}_sort.sam result/STAR/${line1}_L1Aligned.out.sam
    featureCounts -p -O -T 4 -a ${gtfpath}/Mus_musculus.GRCm38.102.chr.gtf -o result/featureCount/${line1}_Counts.txt result/STAR/${line1}_sort.sam

done 3<"SampleName.txt"
