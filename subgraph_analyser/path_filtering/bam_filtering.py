#!/usr/bin/env python3

# Script to filter the results of a read to subgraph analyser path alignment bam file to: 
# 1) identify inverted discordant reads, 2) identify reads with abnormal insert sizes, 
# and 3) identify supplementary, often split, reads to identify insertions that could be due to LGT

import sys
import io
import pysam.samtools
import numpy as np
import collections
import pandas as pd

# Function to filter bamfile for inverted discordant reads
def get_inverted_discordant_reads(bamfile):
    """
    Function to filter bamfile for inverted discordant reads.
    corresponds to samtools flags -f81 -F2048 
    81 being the inverted read flag
    2048 being the secondary alignment flag, which we want to exclude. 
    """
    pysam.samtools.view("-o","inverted_discordant.bam","-f81","-F2048","-@8",bamfile, catch_stdout=False)

def calculate_median_dev(bamfile):
    """
    Function to calculate the median and median absolute deviation of insert sizes from a bam file.
    Adapted from https://github.com/NazifaMoumi/ARGContextProfiler/blob/main/src/context_refinement.py
    """
    #index bamfile if not already indexed
    pysam.samtools.index(bamfile)
    samfile = pysam.AlignmentFile(bamfile, "rb", threads=8)
    all_reads = samfile.fetch()
    size_freq = collections.defaultdict(int)
    for read in all_reads:
        if read.rnext == read.tid and read.is_paired:
            size = abs(read.isize)
            size_freq[size] += 1
            
    all_size = []
    for key, value in size_freq.items():
        # if key (insert size) is 0 skip
        if key == 0:
            continue
        all_size.extend([key] * int(value))
    median_size = np.median(all_size)
    residuals = abs(np.array(all_size) - median_size)
    mad_size = 1.4826 * np.median(residuals)
    return median_size, mad_size

def filter_bam_by_insert_size(bamfile, output_sam_path, median_size, mad_size):
    """
    Filters a BAM file by insert size (TLEN field).

    Args:
        input_bam_path (str): Path to the input BAM file.
        output_bam_path (str): Path to the output filtered BAM file.
        median_size (float): Median insert size.
        mad_size (float): Median absolute deviation of insert sizes.
    """
    try:
        # Define the minimum and maximum insert sizes based on median and MAD
        min_insert_size = median_size - 3 * mad_size
        max_insert_size = median_size + 3 * mad_size
        print(f"Filtering reads with insert sizes less than {min_insert_size} or greater than {max_insert_size}")

        # Open the input BAM file for reading
        infile = pysam.AlignmentFile(bamfile, "rb", threads=8) 

        # Open the output BAM file for writing, using the same header as input
        outfile = pysam.AlignmentFile(output_sam_path, "wb", template=infile)

        # Iterate over each read in the input BAM file
        for read in infile:
            # The 'template_length' attribute (TLEN) represents the insert size
            # It can be negative for the reverse read in a pair
            # Check if the insert size is less than the minimum or greater than the maximum, keep the read
            if read.template_length > max_insert_size:
                outfile.write(read)
            elif read.template_length < (max_insert_size * -1):
                outfile.write(read)
            elif read.template_length < min_insert_size and read.template_length > 0:
                outfile.write(read)
            elif read.template_length > (min_insert_size * -1) and read.template_length < 0:
                outfile.write(read)

        # Close the files
        infile.close()
        outfile.close()
        print(f"Filtered BAM file saved to: {output_sam_path}")

    except FileNotFoundError:
        print(f"Error: Input BAM file not found at {bamfile}")
    except Exception as e:
        print(f"An error occurred: {e}")

# Function to filter bamfile for supplementary reads
def get_supplementary_reads(bamfile):
    """
    Function to filter bamfile for supplementary reads.
    corresponds to samtools flag -f2048 -F1294, supplemental and properly mapped reads
    """
    pysam.samtools.view("-o","supplementary.bam","-f2048","-F1294","-@8",bamfile, catch_stdout=False)


# Function to identify reference paths that meet filtration criteria
def identify_paths_to_remove(inverted_discordant, abnormal_insert_reads, supplementary_reads):
    """
    Function to identify reference paths that meet filtration criteria.
    Filtration criteria:
    inverted reads: inconsistent coverage, specifically <90% or >10% across the path - wrong path missassembly
    abnormal insert size reads: If reads are found at the start (first 10%) or end (last 10%) of the path - syntenic missassembly

    insertion detection criteria:
    abnormal insert size reads: If reads are found in the middle 60% of the path - potential insertion
    supplementary reads: If reads are found in the middle 60% of the path - potential insertion

    output: list of reference paths to remove, list of reference paths with potential insertions
    """
    
    # ## start with inverted discordant reads
    invCoverage_criteria = []
    ## run samtools coverage to get coverage across each path
    inv_coverage_output = pysam.samtools.coverage(inverted_discordant, catch_stdout=True)
    print("calculated coverage for inverted discordant reads")
    # parse coverage output to pandas dataframe
    header = ["rname", "start", "end", "numreads", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq"]
    coverage_df = pd.read_csv(io.StringIO(inv_coverage_output), sep="\t", header=None, names=header, index_col=False, comment="#")
    # separate rname column into subgraph and path at the first underscore
    coverage_df[['subgraph', 'path']] = coverage_df['rname'].str.split('_', n=1, expand=True)

    # number_of_paths_per_subgraph = coverage_df['subgraph'].value_counts().to_dict()
    # multipath_subgraphs = {k: v for k, v in number_of_paths_per_subgraph.items() if v > 1}
    
    #loop through each row in coverage_df to identify those that meet the criteria
    for index, row in coverage_df.iterrows():
        # Add rname to list if coverage is <90% or >10%
        if row['coverage'] < 10 or row['coverage'] > 90:
            invCoverage_criteria.append((row['rname'], 'keep'))
        elif row['coverage'] >= 10 or row['coverage'] <= 90:
            invCoverage_criteria.append((row['rname'], 'discard'))

    ##Abnormal insert size reads
    abnormalInsert_criteria = []
    abnormal_middle_insertion = []
    #index abnormal insert size bam file
    pysam.samtools.index(abnormal_insert_reads)
    # from bamfile get list of reference paths and their lengths
    abnormal_insert_bam = pysam.AlignmentFile(abnormal_insert_reads, "rb", threads=8)
    ref_lengths = {ref: length for ref, length in zip(abnormal_insert_bam.references, abnormal_insert_bam.lengths)}
    abnormal_insert_bam.close()
    
    # loop through ref_lengths
    for ref in sorted(ref_lengths.keys()):
        length = ref_lengths[ref]
        start_abnormal = False
        end_abnormal = False
        header = ["rname", "start", "end", "numreads", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq"]
        #run samtools coverage for the first 5% of the path
        start_5_output = pysam.samtools.coverage("-r", f"{ref}:{0}-{int(length*0.05)}", abnormal_insert_reads, catch_stdout=True)
        # if numreads > 0, add to list with status 'discard'
        start_5_df = pd.read_csv(io.StringIO(start_5_output), sep="\t", header=None, names=header, index_col=False, comment="#")
        if start_5_df['numreads'].iloc[0] > 5:
            start_abnormal = True
        #run samtools coverage for the last 5% of the path
        end_5_output = pysam.samtools.coverage("-r", f"{ref}:{int(length*0.95)}-{length}", abnormal_insert_reads, catch_stdout=True)
        end_5_df = pd.read_csv(io.StringIO(end_5_output), sep="\t", header=None, names=header, index_col=False, comment="#")
        if end_5_df['numreads'].iloc[0] > 5:
            end_abnormal = True
        if start_abnormal and end_abnormal:
            abnormalInsert_criteria.append((ref, 'discard'))
            continue
        #run samtools coverage for the middle 60% of the path
        middle_60_output = pysam.samtools.coverage("-r", f"{ref}:{int(length*0.2)}-{int(length*0.8)}", abnormal_insert_reads, catch_stdout=True)
        middle_60_df = pd.read_csv(io.StringIO(middle_60_output), sep="\t", header=None, names=header, index_col=False, comment="#")
        if middle_60_df['numreads'].iloc[0] > 0:
            abnormal_middle_insertion.append((ref, 'abnormal_reads'))
        else:
            abnormal_middle_insertion.append((ref, 'no_abnormal_reads'))

    ##Supplementary reads
    supplementary_criteria = []
    #index supplementary bam file
    pysam.samtools.index(supplementary_reads)
    # from bamfile get list of reference paths and their lengths
    supplementary_bam = pysam.AlignmentFile(supplementary_reads, "rb", threads=8)
    ref_lengths = {ref: length for ref, length in zip(supplementary_bam.references, supplementary_bam.lengths)}
    supplementary_bam.close()
    # loop through ref_lengths
    for ref in sorted(ref_lengths.keys()):
        length = ref_lengths[ref]
        header = ["rname", "start", "end", "numreads", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq"]
        #run samtools coverage for the middle 60% of the path
        middle_60_output = pysam.samtools.coverage("-r", f"{ref}:{int(length*0.2)}-{int(length*0.8)}", supplementary_reads, catch_stdout=True)
        middle_60_df = pd.read_csv(io.StringIO(middle_60_output), sep="\t", header=None, names=header, index_col=False, comment="#")
        if middle_60_df['numreads'].iloc[0] > 0:
            supplementary_criteria.append((ref, 'supplementary_reads'))
        else:
            supplementary_criteria.append((ref, 'no_supplementary_reads'))
    # combine inverted discordant, abnormal insert size and supplementary criteria into a single dataframe
    inv_df = pd.DataFrame(invCoverage_criteria, columns=['rname', 'inv_status'])
    abnormal_df = pd.DataFrame(abnormalInsert_criteria, columns=['rname', 'abnormal_status'])
    abnormal_middle_df = pd.DataFrame(abnormal_middle_insertion, columns=['rname', 'abnormal_middle_status'])
    supplementary_df = pd.DataFrame(supplementary_criteria, columns=['rname', 'supplementary_status'])
    combined_df = inv_df.merge(abnormal_df, on='rname', how='outer').merge(abnormal_middle_df, on='rname', how='outer').merge(supplementary_df, on='rname', how='outer')
    # fill NaN values in abnormal_status with 'no_abnormal_reads'
    combined_df['abnormal_status'] = combined_df['abnormal_status'].fillna('keep')
    combined_df['abnormal_middle_status'] = combined_df['abnormal_middle_status'].fillna('no_abnormal_reads')
    # fill NaN values in supplementary_status with 'no_supplementary_reads'
    combined_df['supplementary_status'] = combined_df['supplementary_status'].fillna('no_supplementary_reads')
    # write combined_df to csv
    combined_df.to_csv("bam_filtering_summary.csv", index=False)

# __main__
if __name__ == "__main__":
    bamfile=sys.argv[1]
    get_inverted_discordant_reads(bamfile)
    print("Inverted discordant reads extracted to inverted_discordant.bam")
    median_size, mad_size = calculate_median_dev(bamfile)
    print(f"Median insert size: {median_size}, Median absolute deviation: {mad_size}")
    filter_bam_by_insert_size(bamfile, "abnormal_insert.bam", median_size, mad_size)
    print("Abnormal insert size reads extracted to abnormal_insert.bam")
    get_supplementary_reads(bamfile)
    print("Supplementary reads extracted to supplementary.bam")
    identify_paths_to_remove("inverted_discordant.bam", "abnormal_insert.bam", "supplementary.bam")