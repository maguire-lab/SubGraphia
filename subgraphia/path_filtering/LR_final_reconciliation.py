#!/usr/bin/env python3

# Script to take the path metadata from all taxonomically congruent paths and the LR bam filtering summary and reconcile them 
# output: final metadata file, amr gene summary file

import pandas as pd
import sys

def reconcile(path_metadata_file, bam_summary_file, strictness):
    """
    Function to reconcile path metadata with bam filtering summary based on strictness criteria.
    
    :param path_metadata_file: all metadata form path_walk.py
    :param bam_summary_file: LR bam filtering summary from LR_bam_filtering.py
    :param strictness: "strict" or "relaxed" criteria for path retention
    """

    path_metadata = pd.read_csv(path_metadata_file, sep="\t", header=None)
    bam_summary = pd.read_csv(bam_summary_file, sep=",")

    # add metadata header
    path_metadata.columns = ["path_id","path_length","tax_id_lineage","tax_ID_names","LCA_rank"]
    # change rname to path_id in bam summary
    bam_summary.columns = ["path_id", "cov_criterion", "coverage"]
    # join bam summary to path metadata on path_id
    merged = pd.merge(path_metadata, bam_summary, on="path_id", how="right")

    #separate path_id into ARO and path_number on first underscore
    merged[['ARO', 'path_number']] = merged['path_id'].str.split('_', n=1, expand=True)

    # separate path_ID into ARO and path_number on first underscore for path_metadata df
    path_metadata[['ARO', 'path_number']] = path_metadata['path_id'].str.split('_', n=1, expand=True)

    # bin paths into strict or relaxed categories.
    # for LR: strict: only keep paths that meet read alignment criteria (see bam_filtering.py)
    # relaxed: Always keep at least one path per subgraph, even if it does not meet read alignment criteria, take highest coverage path
    #create a dictionary to hold the final status of each path
    final_status = {}
    #create a dictionary counting the number of paths per ARO in path_metadata df, constitutes all the paths that passed taxonomic congruence per ARO
    number_of_paths_per_ARO = path_metadata['ARO'].value_counts().to_dict()
    #Assign final status
    for index, row in merged.iterrows():
        if row['cov_criterion'] == "keep":
            final_status[row['path_id']] = "strict"
        else:
            final_status[row['path_id']] = "discard"
    # find AROs that have no strict paths
    strict_AROs = [key.split('_')[0] for key, value in final_status.items() if value == "strict"]
    all_AROs = list(number_of_paths_per_ARO.keys())
    no_strict_AROs = list(set(all_AROs) - set(strict_AROs))
    #Assign relaxed status to highest coverage path for AROs with no strict paths
    for aro in no_strict_AROs:
        aro_paths = merged[merged['ARO'] == aro]
        # Assign the highest coverage path relaxed status
        highest_cov_path = aro_paths.loc[aro_paths['coverage'].idxmax()]
        final_status[highest_cov_path['path_id']] = "relaxed"
    # create final output dataframe
    final_output = merged.copy()
    final_output['final_status'] = final_output['path_id'].map(final_status)

    # make amr gene summary df
    # find number of paths per ARO after alignment clustering, this is number of paths per ARO in bam summary
    # separate path_id into ARO and path_number on first underscore
    bam_summary[['ARO', 'path_number']] = bam_summary['path_id'].str.split('_', n=1, expand=True)
    number_of_paths_per_ARO_after_clustering = bam_summary['ARO'].value_counts().to_dict()
    
    #find number of paths meeting at least relaxed criteria per ARO
    number_of_relaxed_paths_per_ARO = final_output[final_output['final_status'] != "discard"]['ARO'].value_counts().to_dict()
    # find number of strict paths per ARO
    number_of_strict_paths_per_ARO = final_output[final_output['final_status'] == "strict"]['ARO'].value_counts().to_dict()

    # create amr gene summary dataframe of number_of_paths_per_ARO, number_of_paths_per_ARO_after_clustering, number_of_relaxed_paths_per_ARO
    amr_summary = pd.DataFrame(columns=["ARO", "num_paths_tax_congruent", "num_paths_after_clustering", "num_relaxed_paths", "num_strict_paths"])
    amr_summary['ARO'] = list(number_of_paths_per_ARO.keys())
    amr_summary['num_paths_tax_congruent'] = amr_summary['ARO'].map(number_of_paths_per_ARO)
    amr_summary['num_paths_after_clustering'] = amr_summary['ARO'].map(number_of_paths_per_ARO_after_clustering).fillna(0).astype(int)
    amr_summary['num_relaxed_paths'] = amr_summary['ARO'].map(number_of_relaxed_paths_per_ARO).fillna(0).astype(int)
    amr_summary['num_strict_paths'] = amr_summary['ARO'].map(number_of_strict_paths_per_ARO).fillna(0).astype(int)

    # save amr gene summary file
    amr_summary.to_csv("AMR_genes_summary.tsv", sep="\t", index=False)


    # filter final output based on strictness argument
    if strictness == "strict":
        final_output = final_output[final_output['final_status'] == "strict"]
    elif strictness == "relaxed":
        final_output = final_output[final_output['final_status'].isin(["strict", "relaxed"])]
    else:
        print("Error: strictness argument must be 'strict' or 'relaxed'")
        sys.exit(1)    
    # save final metadata file with path_id	path_length	final_status	tax_id_lineage	tax_ID_names	LCA_rank columns
    final_output = final_output[["path_id", "path_length", "final_status", "tax_id_lineage", "tax_ID_names", "LCA_rank"]]
    final_output.to_csv("final_metadata.tsv", sep="\t", index=False)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 LR_final_reconciliation.py <path_metadata.tsv> <bam_summary.tsv> <strictness>")
        sys.exit(1)
    path_metadata_file = sys.argv[1]
    bam_summary_file = sys.argv[2]
    strictness = sys.argv[3]
    reconcile(path_metadata_file, bam_summary_file, strictness)
