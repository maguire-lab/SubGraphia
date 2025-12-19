#!/usr/bin/env python3

# Script to filter LR bam files based on coverage of the entire path
# Rational: misassembled paths will have gaps of zero coverage, LR graphs are not de Bruijn graphs, edges are not the results of kmer overlaps

import sys
import io
import pysam.samtools
import numpy as np
import collections
import pandas as pd


def filter_clipped_reads(bamfile):
    """
    Function to filter out clipped reads from a LR alignment BAM file. Max allowable clip is 50bp.
    
    :param bamfile: BAM file to be filtered
    :return: Filtered BAM file
    """
    max_clip=50
    pysam.samtools.index(bamfile)
    bam_in = pysam.AlignmentFile(bamfile, "rb")
    bam_out = pysam.AlignmentFile(bamfile.replace(".bam", "_filtered.bam"), "wb", template=bam_in)

    for read in bam_in.fetch():
        cigartuples = read.cigartuples
        if cigartuples is None:
            continue
        # Check for soft or hard clipping at the start or end of the read
        if (cigartuples[0][0] in [4, 5] and cigartuples[0][1] > max_clip) or \
           (cigartuples[-1][0] in [4, 5] and cigartuples[-1][1] > max_clip):
            continue  # Skip clipped reads
        bam_out.write(read)

    bam_in.close()
    bam_out.close()

def identify_paths_to_remove(filtered_bamfile):
    """
    Function to identify paths gaps in coverage in the middle 80% of the path.
    :param filtered_bamfile: Filtered BAM file form filter_clipped_reads
    :return: bam_filtering_summary.csv describing paths to remove and keep 
    """
    header = ["rname", "start", "end", "numreads", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq"]
    # calculate coverage of the middle 80% of each path
    pysam.samtools.index(filtered_bamfile)
    filtered_bam = pysam.AlignmentFile(filtered_bamfile, "rb")
    ref_lengths = {ref: length for ref, length in zip(filtered_bam.references, filtered_bam.lengths)}
    path_coverage = {}
    coverage_df = pd.DataFrame()
    for ref in filtered_bam.references:
        length = ref_lengths[ref]
        start = int(length * 0.1)
        end = int(length * 0.9)
        coverage = pysam.samtools.coverage("-r", f"{ref}:{start}-{end}", filtered_bamfile, catch_stdout=True)
        coverage_df = pd.concat([coverage_df, pd.read_csv(io.StringIO(coverage), sep="\t", header=None, names=header, index_col=False, comment="#")], ignore_index=True)
    filtered_bam.close()
    # Identify paths with less than 98% coverafe in coverage_df
    removal_criteria = []
    # assign keep or discard based on coverage
    for index, row in coverage_df.iterrows():
        if row['coverage'] < 98.0:
            removal_criteria.append((row['rname'], 'discard', row['coverage']))
        else:
            removal_criteria.append((row['rname'], 'keep', row['coverage']))
    removal_df = pd.DataFrame(removal_criteria, columns=['rname', 'coverage_criterion', 'coverage'])
    # save to csv
    removal_df.to_csv("bam_filtering_summary.csv", index=False)

# __main__
if __name__ == "__main__":
    bamfile=sys.argv[1]
    filter_clipped_reads(bamfile)
    identify_paths_to_remove(bamfile.replace(".bam", "_filtered.bam"))