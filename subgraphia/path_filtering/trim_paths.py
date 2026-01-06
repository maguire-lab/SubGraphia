#!/usr/bin/env python3

# Script to find instances where multiple paths assign to the same species, and trimming them to a core aligned region if one exists
import sys
import pandas as pd
from Bio import SeqIO

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

def main(minimap_output_file, metadata_file, amr_gene_summary_file):
    # parse the minimap output
    df = parse_minimap_output(minimap_output_file)
    #filter self alignments
    df = df[df['qname'] != df['tname']]
    # calculate % identity and % query coverage
    mm_out = calculate_id_qc(df)
    # read in metadata file
    metadata_df = pd.read_csv(metadata_file, sep='\t')
    # read in amr gene summary file
    amr_df = pd.read_csv(amr_gene_summary_file, sep='\t')

    # create a dictionary mapping path names to species from the metadata file
    path_to_species = dict(zip(metadata_df['path_id'], metadata_df['tax_ID_names']))
  
    # add qname species and tname species columns to the minimap output dataframe
    mm_out['qname_species'] = mm_out['qname'].map(path_to_species)
    mm_out['tname_species'] = mm_out['tname'].map(path_to_species)

    # add an ARO column for the qname and tname, representing the characters before the first underscore in the path name
    mm_out['qname_ARO'] = mm_out['qname'].str.split('_').str[0]
    mm_out['tname_ARO'] = mm_out['tname'].str.split('_').str[0]

    # filter minimap output to only include alignments where qname and tname have the same species and same ARO
    filtered_mm_out = mm_out[(mm_out['qname_species'] == mm_out['tname_species']) & (mm_out['qname_ARO'] == mm_out['tname_ARO'])]

    # create a set of AROs 
    aro_set = set(filtered_mm_out['qname_ARO'].unique())
    #create empty dataframe to hold paths to trim
    paths_to_trim = {}
    # iterate through the AROs deteministically
    for aro in sorted(aro_set):
        # get the paths associated with this ARO
        aro_df = filtered_mm_out[(filtered_mm_out['qname_ARO'] == aro)]
        # create a species set for this ARO
        species_set = set(aro_df['qname_species'].unique())
        for species in sorted(species_set):
            # filter aro_df to only include this species
            aro_df = aro_df[aro_df['qname_species'] == species]
            # find the max qstart, save the path_ID and qstart
            max_qstart_row = aro_df.loc[aro_df['qstart'].idxmax()]
            max_qstart = max_qstart_row['qstart']
            max_qstart_path = max_qstart_row['qname']
            max_qstart_species = max_qstart_row['qname_species']
            # in rows with same max_qstart_path, find the min qend
            same_start_df = aro_df[aro_df['qname'] == max_qstart_path]
            min_qend_row = same_start_df.loc[same_start_df['qend'].idxmin()]
            min_qend = min_qend_row['qend']
            # save aro, max_qstart_species, path_ID, max_qstart, min_qend to paths_to_trim dataframe
            paths_to_trim[aro] = (max_qstart_species, max_qstart_path, int(max_qstart), int(min_qend))

    # create an ARO column in metadata_df
    metadata_df['ARO'] = metadata_df['path_id'].str.split('_').str[0]

    updated_metadata = []
    aros_trimmed = []
    # loop through metadata_df
    for index, row in metadata_df.iterrows():
        # if path_id ARO is in paths_to_trim, trim the sequence in fastas and write to new fasta file
        aro = row['ARO']
        path_id = row['path_id']
        if aro in paths_to_trim and path_id == paths_to_trim[aro][1] and aro not in aros_trimmed:
            # find the sequence by path_id in fastas
            seq_record = None
            for record in SeqIO.parse(fastas, "fasta"):
                if record.id == path_id:
                    seq_record = record
                    break
            if seq_record is None:
                break
            # trim the sequence
            species, trim_path_id, start, end = paths_to_trim[aro]
            trimmed_seq = seq_record.seq[start:end]
            # write the trimmed sequence to a new fasta file
            trimmed_record = SeqIO.SeqRecord(trimmed_seq, id=path_id+"_trimmed", description="") 
            SeqIO.write(trimmed_record, f"{path_id}_trimmed.fasta", "fasta")
            # update the metadata row with new path_id and length, update final_status to "very_strict"
            row['path_id'] = path_id+"_trimmed"
            row['path_length'] = len(trimmed_seq)
            row['final_status'] = "very_strict"
            updated_metadata.append(row)
            aros_trimmed.append(aro)
        elif aro in paths_to_trim and path_id != paths_to_trim[aro][1]:
            continue
        elif aro in aros_trimmed:
            continue
        else: # aro not in paths_to_trim
            #write out fasta file as is
            seq_record = None
            for record in SeqIO.parse(fastas, "fasta"):
                if record.id == path_id:
                    seq_record = record
                    break
            if seq_record is None:
                break
            SeqIO.write(seq_record, f"{path_id}.fasta", "fasta")
            row['final_status'] = "very_strict"
            updated_metadata.append(row)

    #count paths per ARO in updated_metadata
    aro_path_count = pd.DataFrame(updated_metadata)['ARO'].value_counts().to_dict()
    # Add a 'num_very_strict_paths' column to the amr_df
    amr_df['num_very_strict_paths'] = amr_df['ARO'].map(aro_path_count).fillna(0).astype(int)

    return updated_metadata, amr_df
            





if __name__ == "__main__":
    # specify the input file path
    fastas = sys.argv[1]
    minimap_output_file = sys.argv[2]
    metadata_file = sys.argv[3]
    amr_gene_summary_file = sys.argv[4]

    updated_metadata, amr_df = main(minimap_output_file, metadata_file, amr_gene_summary_file)
    # write updated metadata to a new file
    updated_metadata_df = pd.DataFrame(updated_metadata)
    updated_metadata_df.to_csv("metadata.tsv", sep='\t', index=False)
    # write updated amr gene summary to a new file
    amr_df.to_csv("AMR_genes_summary.tsv", sep='\t', index=False)