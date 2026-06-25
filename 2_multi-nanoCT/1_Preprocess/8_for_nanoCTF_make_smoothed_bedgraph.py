
# ==============================================================================
# Aim:
# Generate smoothed bedGraph tracks from BAM files by counting fragment-end positions
# around paired-end mapped reads.
#
# Process:
#   1. Search each sample folder for:
#        ATACvalidCells_MappedReads_rm_dups_q30_sort.bam
#   2. For each chromosome, keep only properly paired, primary alignments.
#   3. Use one read per fragment to avoid double counting.
#   4. Add signal around both fragment ends using a +/- half_window smoothing window.
#   5. Write non-zero positions as bedGraph:
#        <bam>_smoothed_window<half_window>.bedgraph
#
# Inputs:
#   - <sample>/ATACvalidCells_MappedReads_rm_dups_q30_sort.bam
#
# Outputs:
#   - <sample>/ATACvalidCells_MappedReads_rm_dups_q30_sort_smoothed_window200.bedgraph
#
# Notes:
#   - Edit half_window to change the smoothing window size.
#   - Run this script from the directory containing sample folders.
#   - BAM files must be sorted and indexed before running.
# ==============================================================================

import pysam
import numpy as np
import os

half_window = 200

# make bedgraph
def process_bam(bam_path):
    print(f"Processing: {bam_path}")
    bam = pysam.AlignmentFile(bam_path, "rb")

    for chrom in bam.references:
        chrom_len = bam.get_reference_length(chrom)
        Cuts = np.zeros(chrom_len, dtype='int32')
        pos_to_add = []

        for read in bam.fetch(chrom):
            # Skip unmapped reads or unmapped mates
            if read.is_unmapped or read.mate_is_unmapped:
                continue
            # Skip improper pairs, secondary or supplementary alignments
            if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
                continue
            # Process only reads with positive template length (avoid double counting)
            if read.template_length <= 0:
                continue
            start = read.reference_start
            end = read.next_reference_start + abs(read.template_length)
            if start >= end or end > chrom_len:
                continue
            # Collect positions to add coverage around start and end with smoothing window
            pos_to_add.extend(range(max(0, start - half_window), min(start + half_window + 1, chrom_len)))
            pos_to_add.extend(range(max(0, end - half_window), min(end + half_window + 1, chrom_len)))

        # Add coverage counts in bulk for performance
        np.add.at(Cuts, pos_to_add, 1)

        # Extract non-zero coverage positions
        nonzero_pos = np.nonzero(Cuts)[0]
        output_lines = [f"{chrom}\t{i}\t{i+1}\t{Cuts[i]}" for i in nonzero_pos]

        output_path = bam_path.replace(".bam", f"_smoothed_window{half_window}.bedgraph")
        # Write coverage to bedgraph file
        with open(output_path, "a") as f_out:
            f_out.write("\n".join(output_lines) + "\n")

    bam.close()
    print(f"Finished: {output_path}")


if __name__ == "__main__":
    # List directories in current working directory
    folders = [f for f in os.listdir() if os.path.isdir(f)]
    print(f"\nFound folders: {folders}")

    for folder in folders:
        bam_files = [
            # os.path.join(folder, "Mapped_rm_dups_q30_sort.bam"),
            os.path.join(folder, "ATACvalidCells_MappedReads_rm_dups_q30_sort.bam")
        ]

        for bam_path in bam_files:
            if os.path.exists(bam_path):
                process_bam(bam_path)
            else:
                print(f"File not found: {bam_path}")
