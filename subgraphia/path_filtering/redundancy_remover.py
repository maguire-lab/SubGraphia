#!/usr/bin/env python3

# Script to parse the results of an all vs all minimap alignment of fasta files (ava_minimap.sh output) and identify fastas that are highly redundant to longer fastas
import sys
import pandas as pd

# define a function to parse the minimap output
def parse_minimap_output(file_path):
    # read the minimap output into a pandas DataFrame
    col_names=['qname', 'qlen', 'qstart', 'qend', 'strand', 'tname', 'tlen', 'tstart', 'tend', 'nmatch', 'alnlen']
    df = pd.read_csv(file_path, sep='\t', header=None)
    #keep only the first 11 columns
    df = df.iloc[:, :11]
    df.columns = col_names
    return df

# define a function to calculate % identity and % query coverage from the minimap output
def calculate_id_qc(df):
    df['id'] = df['nmatch'] / df['alnlen'] * 100
    df['qc'] = (df['qend'] - df['qstart']) / df['qlen'] * 100
    return df

# define a function to remove redundant sequences
def remove_redundant_sequences(df, id_threshold, qc_threshold):
    # sort sequences longest to shortest
    seq_lens = df.groupby('qname')['qlen'].max()
    seq_lens = seq_lens.sort_values(ascending=False)
    #convert seq_lens to a list of just the names
    seq_names = seq_lens.index.tolist()

    redundant_seqs= []
    representative_seqs = []

    #filter away self alignments
    df = df[df['qname'] != df['tname']]
    
    # iterate through the sequences, bin into redundant and representative sequences
    for seq in seq_names:
        if seq in redundant_seqs:
            continue
        # get the sequences that are redundant to the current sequence
        redundant = df[(df['tname'] == seq) & (df['id'] >= id_threshold) & (df['qc'] >= qc_threshold)]
        #append all qnames in redundant to the redundant_seqs list
        redundant_seqs.extend(redundant['qname'].tolist())
        #Add the current sequence to the representative sequences
        representative_seqs.append(seq)

    return representative_seqs

if __name__ == "__main__":
    # specify the input file path
    minimap_output_file = sys.argv[1]
    id_threshold = int(sys.argv[2])  # identity threshold
    qc_threshold = int(sys.argv[3])  # query coverage threshold

    # parse the minimap output
    df = parse_minimap_output(minimap_output_file)

    # calculate % identity and % query coverage
    df = calculate_id_qc(df)

    # remove redundant sequences
    representative_seqs = remove_redundant_sequences(df, id_threshold, qc_threshold)
    
    # write the representative sequences to a file with .fasta appended to each item in the list
    with open('representative_sequences.tsv', 'w') as f:
        for seq in representative_seqs:
            f.write(f"{seq}\n")
