#!/usr/bin/env python3

# Script to take the path metadata from all taxonomically congruent paths and the bam filtering summary and reconcile them 
# output: final metadata file, amr gene summary file
# usage: python3 final_reconciliation.py <path_metadata.tsv> <bam_summary.tsv> <strictness>
# strict definition: only keep paths that meet read alignment criteria (see bam_filtering.py)
# relaxed definition: Always keep at least one path per subgraph, even if it does not meet read alignment criteria

import pandas as pd
import sys

def reconcile(path_metadata_file, bam_summary_file, strictness):
    # read in the path metadata and bam summary files
    path_metadata = pd.read_csv(path_metadata_file, sep="\t", header=None)
    bam_summary = pd.read_csv(bam_summary_file, sep=",")

    # add metadata header
    path_metadata.columns = ["path_id", "length", "taxid_lineage", "lineage", "coverage", "coverage_stdev"]
    # change rname to path_id in bam summary
    bam_summary.columns = ["path_id", "inv_status", "abnormal_status", "abnormal_middle_status", "supplementary_status", "clipped_status"]

    # join bam summary to path metadata on path_id
    merged = pd.merge(path_metadata, bam_summary, on="path_id", how="right")

    #separate path_id into ARO and path_number on first underscore
    merged[['ARO', 'path_number']] = merged['path_id'].str.split('_', n=1, expand=True)

    # separate path_ID into ARO and path_number on first underscore for path_metadata df
    path_metadata[['ARO', 'path_number']] = path_metadata['path_id'].str.split('_', n=1, expand=True)

    # bin paths into strict or relaxed categories. 
    ## criteria for strict: inv_status == "keep" and abnormal_status == "keep"
    ## criteria for relaxed: if no paths of a given ARO meet strict criteria, keep the path that satisfies the inv_status == "keep" criteria, ifelse the abnormal_status == "keep" criteria, ifelse keep the longest path
    ## paths not strict or relaxed are labeled "discard"
    
    #create a dictionary to hold the final status of each path
    final_status = {}

    #create a dictionary counting the number of paths per ARO in path_metadata df, constitutes all the paths that passed taxonomic congruence per ARO
    number_of_paths_per_ARO = path_metadata['ARO'].value_counts().to_dict()

    #Assign final status
    for index, row in merged.iterrows():
        if row['inv_status'] == "keep" and row['abnormal_status'] == "keep" and row['clipped_status'] == "keep":
            final_status[row['path_id']] = "strict"
        elif row['inv_status'] == "keep" or row['abnormal_status'] == "keep" or row['clipped_status'] == "keep":
            final_status[row['path_id']] = "relaxed"
        else: # inv_status, abnormal_status, and clipped_status are "discard"
            final_status[row['path_id']] = "discard"

    # find AROs that have no strict or relaxed paths
    strict_or_relaxed_AROs = [key.split('_')[0] for key, value in final_status.items() if value in ["strict", "relaxed"]]
    all_AROs = list(number_of_paths_per_ARO.keys())
    no_strict_or_relaxed_AROs = list(set(all_AROs) - set(strict_or_relaxed_AROs))

    #Assign relaxed status to longest path for AROs with no strict or relaxed paths
    for aro in no_strict_or_relaxed_AROs:
        aro_paths = merged[merged['ARO'] == aro]
        # Assign the longest path relaxed status
        longest_path = aro_paths.loc[aro_paths['length'].idxmax()]
        final_status[longest_path['path_id']] = "relaxed"

    # add final status to merged dataframe
    merged['final_status'] = merged['path_id'].map(final_status)

    # summarise the final_status counts and print
    print(merged['final_status'].value_counts())

    # filter merged dataframe to only include paths that meet the specified strictness criteria
    if strictness == "strict":
        final_df = merged[merged['final_status'] == "strict"]
    elif strictness == "relaxed":
        final_df = merged[merged['final_status'].isin(["strict", "relaxed"])]
    else:
        print("Error: strictness must be 'strict' or 'relaxed'")
        sys.exit(1)

    # make final_df nicer
    final_df = final_df[["path_id", 
                         "length", 
                         "final_status", 
                        #  "abnormal_middle_status", # Could potentially include this column again in future, could be used for LGT detection
                        #  "supplementary_status", 
                         "coverage", 
                         "coverage_stdev", 
                         "taxid_lineage", 
                         "lineage"]]

    # make amr gene summary df
    # find number of paths per ARO after alignment clustering, this is number of paths per ARO in bam summary
    # separate path_id into ARO and path_number on first underscore
    bam_summary[['ARO', 'path_number']] = bam_summary['path_id'].str.split('_', n=1, expand=True)
    number_of_paths_per_ARO_after_clustering = bam_summary['ARO'].value_counts().to_dict()
    
    #find number of paths meeting at least relaxed criteria per ARO
    number_of_relaxed_paths_per_ARO = merged[merged['final_status'] != "discard"]['ARO'].value_counts().to_dict()
    # find number of strict paths per ARO
    number_of_strict_paths_per_ARO = merged[merged['final_status'] == "strict"]['ARO'].value_counts().to_dict()

    # create amr gene summary dataframe of number_of_paths_per_ARO, number_of_paths_per_ARO_after_clustering, number_of_relaxed_paths_per_ARO
    amr_summary = pd.DataFrame(columns=["ARO", "num_paths_tax_congruent", "num_paths_after_clustering", "num_relaxed_paths", "num_strict_paths"])
    amr_summary['ARO'] = list(number_of_paths_per_ARO.keys())
    amr_summary['num_paths_tax_congruent'] = amr_summary['ARO'].map(number_of_paths_per_ARO)
    amr_summary['num_paths_after_clustering'] = amr_summary['ARO'].map(number_of_paths_per_ARO_after_clustering).fillna(0).astype(int)
    amr_summary['num_relaxed_paths'] = amr_summary['ARO'].map(number_of_relaxed_paths_per_ARO).fillna(0).astype(int)
    amr_summary['num_strict_paths'] = amr_summary['ARO'].map(number_of_strict_paths_per_ARO).fillna(0).astype(int)

    # sort final_df by path_id
    final_df = final_df.sort_values(by="path_id")

    # write final_df to tsv
    final_df.to_csv("final_metadata.tsv", sep="\t", index=False)

    # Sort amr_summary by ARO
    amr_summary = amr_summary.sort_values(by="ARO")

    # write amr_summary to tsv
    amr_summary.to_csv("AMR_genes_summary.tsv", sep="\t", index=False)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 final_reconciliation.py <path_metadata.tsv> <bam_summary.tsv> <strictness>")
        sys.exit(1)
    path_metadata_file = sys.argv[1]
    bam_summary_file = sys.argv[2]
    strictness = sys.argv[3]
    reconcile(path_metadata_file, bam_summary_file, strictness)
