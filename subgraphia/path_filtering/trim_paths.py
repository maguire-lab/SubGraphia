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
    # create empty final multifasta list to hold the final sequences that can be written out after this function is called
    final_multifasta = []
    # parse the minimap output
    df = parse_minimap_output(minimap_output_file)
    #filter self alignments
    df = df[df['qname'] != df['tname']]
    # calculate % identity and % query coverage
    mm_out = calculate_id_qc(df)
    # filter away alignments with less than 90% identity
    mm_out = mm_out[mm_out['id'] >= 90]
    # filter away alignments with less than 15% query coverage
    mm_out = mm_out[mm_out['qc'] >= 15]

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
            aro_df_filt = aro_df[aro_df['qname_species'] == species]
            # find the max qstart, save the path_ID and qstart
            max_qstart_row = aro_df_filt.loc[aro_df_filt['qstart'].idxmax()]
            max_qstart = max_qstart_row['qstart']
            max_qstart_path = max_qstart_row['qname']
            max_qstart_species = max_qstart_row['qname_species']
            # in rows with same max_qstart_path, find the min qend
            same_start_df = aro_df_filt[aro_df_filt['qname'] == max_qstart_path]
            min_qend_row = same_start_df.loc[same_start_df['qend'].idxmin()]
            min_qend = min_qend_row['qend']
            if min_qend <= max_qstart:
                # find second min qend that is greater than max_qstart, can't be equal to max_qstart because that would mean the path is being trimmed to 0 length
                min_qend_row = same_start_df[same_start_df['qend'] > max_qstart].loc[same_start_df[same_start_df['qend'] > max_qstart]['qend'].idxmin()]
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
        species = row['tax_ID_names']
        if aro in paths_to_trim and species == paths_to_trim[aro][0] and path_id == paths_to_trim[aro][1] and aro not in aros_trimmed:
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
            # Add the trimmed sequence to the final multifasta list as a SeqRecord object with id = path_id+"_trimmed" and description = ""
            trimmed_seq_record = SeqIO.SeqRecord(trimmed_seq, id=path_id+"_trimmed", description="")
            final_multifasta.append(trimmed_seq_record)
            # update the metadata row with new path_id and length, update final_status to "very_strict"
            row['path_id'] = path_id+"_trimmed"
            row['path_length'] = len(trimmed_seq)
            row['final_status'] = "very_strict"
            updated_metadata.append(row)
            aros_trimmed.append(aro)
        elif aro in paths_to_trim and species == paths_to_trim[aro][0] and path_id != paths_to_trim[aro][1]:
            continue
        elif aro in aros_trimmed and species == paths_to_trim[aro][0]:
            continue
        else: # aro not in paths_to_trim
            # Add the original sequence to the final multifasta list as a SeqRecord object with id = path_id and description = ""
            seq_record = None
            for record in SeqIO.parse(fastas, "fasta"):
                if record.id == path_id:
                    seq_record = record
                    break
            if seq_record is None:
                break
            final_multifasta.append(seq_record)
            row['final_status'] = "very_strict"
            updated_metadata.append(row)

    #count paths per ARO in updated_metadata
    aro_path_count = pd.DataFrame(updated_metadata)['ARO'].value_counts().to_dict()
    # Add a 'num_very_strict_paths' column to the amr_df
    amr_df['num_very_strict_paths'] = amr_df['ARO'].map(aro_path_count).fillna(0).astype(int)

    return final_multifasta, updated_metadata, amr_df

def pathSeq_v_geneSeq(final_multifasta, updated_metadata, amr_df, gene_seqs):
    """
    Function to compare the trimmed path sequences to the reference gene sequences, the following rules apply:
    1) AROs with 0 paths are assigned the gene sequence
    2) AROs with >= 1 path shorter than the gene sequence are assigned the gene sequence, other paths removed
    3) AROs with some shorter paths and some longer paths than the gene sequence have the shorter paths removed and the longer paths kept, gene sequence not added
    4) AROs with all paths longer than the gene sequence, nothing changed. 

    returns:
    multifasta: modified as above
    metadata: modified as above, new path_IDs, new path lengths, NAs for tax_ID_names, tax_ID_lineage, and LCA_Rank
    AMR_summary: modified as above, updated path counts per ARO
    """

    # parse gene_seqs fasta file into a dictionary mapping ARO to sequence 
    gene_seq_dict = {}
    for record in SeqIO.parse(gene_seqs, "fasta"):
        gene_seq_dict[record.id] = record.seq

    # create a dictionary mapping ARO to a list of path sequence lengths from updated_metadata
    aro_to_path_lengths = {}
    for row in updated_metadata:
        aro = row['ARO']
        path_length = row['path_length']
        # using AROs as keys, append all path lengths to a list for each ARO
        if aro in aro_to_path_lengths:
            aro_to_path_lengths[aro].append(path_length)
        else:            
            aro_to_path_lengths[aro] = [path_length]

    # create empty lists to hold the final multifasta, updated metadata, and updated amr_df
    multifasta = []
    metadata = []
    amr_summary = []

    # loop through amr_df, determine number of paths per ARO, get lengths of paths from updated_metadata, get gene sequence lengths from gene_seqs, apply rules above
    for index, row in amr_df.iterrows():
        aro = row['ARO']
        num_paths = row['num_very_strict_paths']
        # if num_paths == 0
        if num_paths == 0:
            # add the gene sequence to the multifasta with id = aro+"_gene_seq" and description = ""
            gene_seq_record = SeqIO.SeqRecord(gene_seq_dict[aro], id=aro+"_gene_seq", description="")
            multifasta.append(gene_seq_record)
            # add a row to metadata with path_id = aro+"_gene_seq", path_length = length of gene sequence, tax_ID_names = NA, tax_ID_lineage = NA, LCA_Rank = NA, ARO = aro, final_status = "very_strict"
            metadata.append({'path_id': aro+"_gene_seq", 'path_length': len(gene_seq_dict[aro]), 'tax_ID_names': "NA", 'tax_id_lineage': "NA", 'LCA_rank': "NA", 'ARO': aro, 'final_status': "very_strict"})
            # add a row to amr_summary with ARO = aro, num_very_strict_paths = 1
            amr_summary.append({'ARO': aro, 'num_very_strict_paths': 1})
        elif num_paths >= 1:
            path_lengths = aro_to_path_lengths[aro]
            gene_seq_length = len(gene_seq_dict[aro])
            # if all paths are shorter than the gene sequence, add the gene sequence to the multifasta with id = aro+"_gene_seq" and description = "", do not add any of the existing paths to the updated files
            if all(path_length < gene_seq_length for path_length in path_lengths):
                gene_seq_record = SeqIO.SeqRecord(gene_seq_dict[aro], id=aro+"_gene_seq", description="")
                multifasta.append(gene_seq_record)
                metadata.append({'path_id': aro+"_gene_seq", 'path_length': len(gene_seq_dict[aro]), 'tax_ID_names': "NA", 'tax_id_lineage': "NA", 'LCA_rank': "NA", 'ARO': aro, 'final_status': "very_strict"})
                amr_summary.append({'ARO': aro, 'num_very_strict_paths': 1})
            # if some paths are shorter than the gene sequence and some paths are longer than the gene sequence, add the longer paths to the multifasta and updated metadata, do not add the gene sequence
            elif any(path_length < gene_seq_length for path_length in path_lengths) and any(path_length >= gene_seq_length for path_length in path_lengths):
                for row in updated_metadata:
                    if row['ARO'] == aro and row['path_length'] >= gene_seq_length:
                        # add the path sequence to the multifasta with id = path_id and description = ""
                        seq_record = None
                        for record in final_multifasta:
                            if record.id == row['path_id']:
                                seq_record = record
                                break
                        if seq_record is None:
                            break
                        multifasta.append(seq_record)
                        # add a row to metadata with path_id = path_id, path_length = path_length, tax_ID_names = tax_ID_names, tax_ID_lineage = tax_ID_lineage, LCA_Rank = LCA_Rank, ARO = aro, final_status = final_status from updated_metadata
                        metadata.append({'path_id': row['path_id'], 'path_length': row['path_length'], 'tax_ID_names': row['tax_ID_names'], 'tax_id_lineage': row['tax_id_lineage'], 'LCA_rank': row['LCA_rank'], 'ARO': aro, 'final_status': row['final_status']})
                amr_summary.append({'ARO': aro, 'num_very_strict_paths': sum(1 for path_length in path_lengths if path_length >= gene_seq_length)})
            # if all paths are longer than the gene sequence, add all paths to the multifasta and updated metadata, do not add the gene sequence
            elif all(path_length >= gene_seq_length for path_length in path_lengths):
                for row in updated_metadata:
                    if row['ARO'] == aro:
                        # add the path sequence to the multifasta with id = path_id and description = ""
                        seq_record = None
                        for record in final_multifasta:
                            if record.id == row['path_id']:
                                seq_record = record
                                break
                        if seq_record is None:
                            break
                        multifasta.append(seq_record)
                        # add a row to metadata with path_id = path_id, path_length = path_length, tax_ID_names = tax_ID_names, tax_ID_lineage = tax_ID_lineage, LCA_Rank = LCA_Rank, ARO = aro, final_status = final_status from updated_metadata
                        metadata.append({'path_id': row['path_id'], 'path_length': row['path_length'], 'tax_ID_names': row['tax_ID_names'], 'tax_id_lineage': row['tax_id_lineage'], 'LCA_rank': row['LCA_rank'], 'ARO': aro, 'final_status': row['final_status']})
                amr_summary.append({'ARO': aro, 'num_very_strict_paths': num_paths})

    # create a new amr_df from amr_summary
    amr_summary_df = pd.DataFrame(amr_summary)

    # merge amr_summary_df with original amr_df to get all the original columns back, using ARO as the key
    # drop num_very_strict_paths column from original amr_df
    to_join_amr_df = amr_df.drop(columns=['num_very_strict_paths'])
    merged_amr_df = pd.merge(to_join_amr_df, amr_summary_df, on='ARO', how='left')

    # create a new metadata_df from metadata
    metadata_df = pd.DataFrame(metadata)

    return multifasta, metadata_df, merged_amr_df

if __name__ == "__main__":
    # specify the input file path
    fastas = sys.argv[1]
    minimap_output_file = sys.argv[2]
    metadata_file = sys.argv[3]
    amr_gene_summary_file = sys.argv[4]
    gene_seqs = sys.argv[5]
    final_multifasta, updated_metadata, amr_df = main(minimap_output_file, metadata_file, amr_gene_summary_file)
    multifasta, metadata_df, merged_amr_df = pathSeq_v_geneSeq(final_multifasta, updated_metadata, amr_df, gene_seqs)
    # write each sequence in the final multifasta list to a separate fasta file with the name of the sequence ID
    for seq_record in multifasta:
        SeqIO.write(seq_record, f"{seq_record.id}.fasta", "fasta")
    # write updated metadata to a new file
    metadata_df.to_csv("metadata.tsv", sep='\t', index=False)
    # write updated amr summary to a new file
    merged_amr_df.to_csv("AMR_genes_summary.tsv", sep='\t', index=False)