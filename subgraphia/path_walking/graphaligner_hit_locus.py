#!/usr/bin/env python3

## script to process the output of graphaligner run against CARD to assign only on ehit to each graph location
## USAGE: python3 graphaligner_hit_locus.py <graphaligner_INPUT> <output_file.tsv>



import sys
import pandas as pd

def process_graphaligner_output(graphaligner_output):
    # define column names
    graphaligner_cols = [
        "query", "query_length", "query_start", "query_end", "strand", "path", 
        "path_length", "path_start", "path_end", "number_residue_matches", 
        "alignment_length", "mapping_quality", "snps_and_indels", 
        "alignment_score", "divergence", "perc_identity", "cigar"
    ]

    # read in the graphaligner output (tsv)
    graphaligner_df = pd.read_csv(graphaligner_output, sep="\t", header=None, names=graphaligner_cols)

    # remove id:f: from the per_identity column
    graphaligner_df["perc_identity"] = graphaligner_df["perc_identity"].str.replace("id:f:", "")
    # convert to float
    graphaligner_df["perc_identity"] = graphaligner_df["perc_identity"].astype(float)

    # calculate percent query coverage
    graphaligner_df["perc_query_cov"] = (graphaligner_df["query_end"] - graphaligner_df["query_start"])/graphaligner_df["query_length"]

    # separate out query column 
    graphaligner_df[["query_name", "ARO", "family"]] = graphaligner_df["query"].str.split("|", expand=True)

    # filter away any hits less than 50% query coverage or 50% identity
    graphaligner_df = graphaligner_df[graphaligner_df["perc_query_cov"] >= 0.5]
    graphaligner_df = graphaligner_df[graphaligner_df["perc_identity"] >= 0.5]

    # create a list of unique paths in the graphaligner output
    paths = graphaligner_df["path"].unique()

    ## identify the hit loci 
    graphaligner_df["hit_region"] = 0

    # initialize hit region counter
    counter = 0

    # for each unique segment in the graphaligner results
    for seg in paths:
        # for each instance of that segment
        for idx in graphaligner_df[graphaligner_df["path"] == seg].index:
            # Find that path start where a hit region hasn't been assigned, smallest value on segment
            PS = graphaligner_df[(graphaligner_df["path"] == seg) & (graphaligner_df["hit_region"] == 0)]["path_start"].min()
            
            # find the max alignment length for that PS
            aln_max = graphaligner_df[(graphaligner_df["path"] == seg) & (graphaligner_df["hit_region"] == 0) & (graphaligner_df["path_start"] == PS)]["alignment_length"].max()
            
            # path end calculation
            PE = PS + aln_max
            
            # if variables have successfully been identified
            if not pd.isna(PS) and not pd.isna(PE):
                # add a new hit region number
                counter += 1
                
                # for each row in graphaligner_df, look for hit region coordinates and assign new hit region
                for n in graphaligner_df.index:
                    if (graphaligner_df.at[n, "path"] == seg and 
                        graphaligner_df.at[n, "path_start"] >= PS and 
                        graphaligner_df.at[n, "path_end"] <= PE and 
                        graphaligner_df.at[n, "hit_region"] == 0):
                        graphaligner_df.at[n, "hit_region"] = counter
                
                # remove used coordinates
                del PS
                del PE
    
    # Summarize data for each hit region
    summarized_hits = graphaligner_df.sort_values(
        by=["perc_query_cov", "perc_identity", "alignment_length", "ARO"], 
        ascending=[False, False, False, True]
    ).drop_duplicates(subset=["hit_region"], keep="first")

    return summarized_hits

if __name__ == '__main__':
    graphaligner_output = sys.argv[1]
    summarized_hits = process_graphaligner_output(graphaligner_output)
    # must make ARO column unique. Add _instance_1,2,3 etc as a suffix
    summarized_hits["ARO"] = summarized_hits.groupby("ARO").cumcount().add(1).astype(str).radd(summarized_hits["ARO"] + "inst")
    summarized_hits.to_csv(sys.argv[2], sep="\t", index=False)