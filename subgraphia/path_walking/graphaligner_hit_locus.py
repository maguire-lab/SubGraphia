#!/usr/bin/env python3

## script to process the output of graphaligner run against CARD to assign only on ehit to each graph location
## USAGE: python3 graphaligner_hit_locus.py <graphaligner_INPUT> <output_file.tsv>

import networkx as nx
import gfapy
import sys
import re
import os
import pandas as pd
from Bio.Seq import Seq

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

    # separate out query column "gb|AF028812.1|+|392-887|ARO:3002867|dfrF [Enterococcus faecalis]"
    graphaligner_df[["db","accn","strand","coords","ARO","gene"]] = graphaligner_df["query"].str.split("|", expand=True)
    # extract ARO number only
    graphaligner_df["ARO"] = graphaligner_df["ARO"].str.split(":").str[1]

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

def parse_to_directed_graph(gfa_input):
  gfa_graph = gfapy.Gfa.from_file(gfa_input)
  directed_graph = nx.DiGraph()

  ### create nodes
  for segment in gfa_graph.segments:
    directed_graph.add_node(segment.name, name = segment.name, sequence = segment.sequence, KC = segment.KC, DP = segment.DP)

  #get edges from gfa_graph
  edges = [line for line in gfa_graph.lines if isinstance(line, gfapy.Line) and line.record_type == "L"]

  # get overlap length from edges
  if edges:
    overlap = int(str(edges[0].overlap)[:-1])
  else:
    overlap = 0

  for edge in edges:
    directed_graph.add_edge(edge.from_segment.name, edge.to_segment.name, from_orient = edge.from_orient, to_orient = edge.to_orient)
  
  return directed_graph, overlap

def get_path_sequences(summarized_hits, gfa_input):
    """
    Function to take the paths of the graphaligner output and return the sequences of those paths from the directed graph. 
    Directed graph overlap must be supplied, intended for use with graph that have an overlap of k and a minimum segment length of k+1
    """
    # create tuple with ARO and path
    congruent_paths = list(zip(summarized_hits["ARO"], summarized_hits["path"]))
    # remove > and < beginning and end of path
    congruent_paths = [(ARO, path.strip("<>")) for ARO, path in congruent_paths]
    # split path into list of nodes, replace internal > and < with ","
    congruent_paths = [(ARO, re.split(r'[><]', path)) for ARO, path in congruent_paths]

    #create an empty multifasta
    fasta_content = ""

    # output a fasta file for each path 
    for ARO, path in congruent_paths:
        print("Processing path for ARO: " + ARO + " at path: " + ">" + ">".join(path) + "<")
        # use regular expressions to subset the gfa input to only the segments and links that are in the path
        subset_gfa = ""
        with open(gfa_input, 'r') as gfa_file:
            for line in gfa_file:
                if re.search(r'^[SL]\t(' + '|'.join(path) + r')\b', line):
                    subset_gfa += line
        # if any link lines have an edge that is not in a segment line, add a placeholder empty segment line for that node to the subset gfa
        subset_gfa_lines = subset_gfa.splitlines()
        existing_segments = set()
        for line in subset_gfa_lines:
            if line.startswith("S"):
                seg_name = line.split("\t")[1]
                existing_segments.add(seg_name)
        # Track segments we add as placeholders to avoid duplicates
        added_segments = set()
        for line in subset_gfa_lines:
            if line.startswith("L"):
                from_node = line.split("\t")[1]
                to_node = line.split("\t")[3]
                # Check if from_node is present in any S line or already added
                if from_node not in existing_segments and from_node not in added_segments:
                    subset_gfa += f"S\t{from_node}\t*\n"
                    added_segments.add(from_node)
                # Check if to_node is present in any S line or already added
                if to_node not in existing_segments and to_node not in added_segments:
                    subset_gfa += f"S\t{to_node}\t*\n"
                    added_segments.add(to_node)
        
        # save subset gfa to a file
        with open(gfa_input.split(".")[0] + "_" + ARO + "_subset.gfa", 'w') as subset_gfa_file:
            subset_gfa_file.write(subset_gfa)
        # load the subset gfa into a directed graph
        directed_graph, overlap = parse_to_directed_graph(gfa_input.split(".")[0] + "_" + ARO + "_subset.gfa")
        # if path contains only one node, add the sequence of that node to fasta content
        if len(path) == 1:
            node = path[0]
            node_data = directed_graph.nodes[node]
            fasta_header = ">" + ARO
            fasta_content += fasta_header + "\n" + node_data['sequence'] + "\n"
            continue
        
        for node in path:
            # Establish a directionality for the path, if an edge has the same directionality, no reverse complement is needed, if it has a different directionality, reverse complement is needed
            # start with the first node, get the sequence and orientation
            if node == path[0]:
                #if out edge:
                if directed_graph.get_edge_data(node, path[path.index(node) + 1]) is not None:
                    path_dir = "out"
                    out_edge_data = directed_graph.get_edge_data(node, path[path.index(node) + 1])
                # if in edge:
                elif directed_graph.get_edge_data(path[path.index(node) + 1], node) is not None:
                    path_dir = "in"
                    in_edge_data = directed_graph.get_edge_data(path[path.index(node) + 1], node)
                else:
                    print("No edge found between " + node + " and " + path[path.index(node) + 1])
                    break
                # Get the sequence 
                node_data = directed_graph.nodes[node]
                if path_dir == "out":
                    # get the sequence orientation from out edge data
                    node_orient = out_edge_data['from_orient']
                else:
                    # get the sequence orientation from in edge data
                    node_orient = in_edge_data['to_orient']
                # trim overlap from the right end of the sequence
                if node_orient == "+": #do not reverse complement
                    if overlap >0:
                        trimmed_sequence = node_data['sequence'][:-overlap]
                    else:
                        trimmed_sequence = node_data['sequence']
                    path_sequence = trimmed_sequence
                elif node_orient == "-": #reverse complement
                    RC_sequence = Seq(node_data['sequence']).reverse_complement()
                    if overlap >0:
                        trimmed_sequence = str(RC_sequence)[:-overlap]
                    else:
                        trimmed_sequence = str(RC_sequence)
                    path_sequence = trimmed_sequence

            # skip if at the last node
            if node == path[-1]:
                continue

            # now add the next node in the path for each node
            #check to see if edge matches the directionality of the path
            if directed_graph.get_edge_data(node, path[path.index(node) + 1]) is not None:
                node_dir = "out"
                out_edge_data = directed_graph.get_edge_data(node, path[path.index(node) + 1])
                node_orient = out_edge_data['to_orient']
            elif directed_graph.get_edge_data(path[path.index(node) + 1], node) is not None:
                node_dir = "in"
                in_edge_data = directed_graph.get_edge_data(path[path.index(node) + 1], node)
                node_orient = in_edge_data['from_orient']
            else:
                print("No edge found between " + node + " and " + path[path.index(node) + 1])
                break
            if path_dir == node_dir:
                # directionality matches, reverse complement as normal
                node_data = directed_graph.nodes[path[path.index(node) + 1]]
                if node_orient == "+":
                    if overlap >0:
                        trimmed_sequence = node_data['sequence'][:-overlap]
                    else:
                        trimmed_sequence = node_data['sequence']
                    #if the path direction is "out" add trimmed sequence to the end of the path sequence
                    if path_dir == "out":
                        path_sequence += trimmed_sequence
                    #if the path direction is "in" add trimmed sequence to the start of the path sequence
                    elif path_dir == "in":
                        path_sequence = trimmed_sequence + path_sequence
                elif node_orient == "-":
                    RC_sequence = Seq(node_data['sequence']).reverse_complement()
                    if overlap >0:
                        trimmed_sequence = str(RC_sequence)[:-overlap]
                    else:
                        trimmed_sequence = str(RC_sequence)
                    #if the path direction is "out" add trimmed sequence to the end of the path sequence
                    if path_dir == "out":
                        path_sequence += trimmed_sequence
                    #if the path direction is "in" add trimmed sequence to the start of the path sequence
                    elif path_dir == "in":
                        path_sequence = trimmed_sequence + path_sequence
            elif path_dir != node_dir:
                # directionality does not match, reverse reverse complementing 
                node_data = directed_graph.nodes[path[path.index(node) + 1]]
                if node_orient == "+": 
                    RC_sequence = Seq(node_data['sequence']).reverse_complement()
                    if overlap >0:
                        trimmed_sequence = str(RC_sequence)[:-overlap]
                    else:
                        trimmed_sequence = str(RC_sequence)
                    #if the path direction is "out" add trimmed sequence to the end of the path sequence
                    if path_dir == "out":
                        path_sequence += trimmed_sequence
                    #if the path direction is "in" add trimmed sequence to the start of the path sequence
                    elif path_dir == "in":
                        path_sequence = trimmed_sequence + path_sequence
                elif node_orient == "-":
                    if overlap >0:
                        trimmed_sequence = node_data['sequence'][:-overlap]
                    else:
                        trimmed_sequence = node_data['sequence']
                    #if the path direction is "out" add trimmed sequence to the end of the path sequence
                    if path_dir == "out":
                        path_sequence += trimmed_sequence
                    #if the path direction is "in" add trimmed sequence to the start of the path sequence
                    elif path_dir == "in":
                        path_sequence = trimmed_sequence + path_sequence

        fasta_header = ">" + ARO
        #write the fasta header and sequence to the fasta content
        fasta_content += fasta_header + "\n" + path_sequence + "\n"
    
    # write the fasta content to a file
    with open(gfa_input.split(".")[0] + "_gene_sequences.fasta", 'w') as fasta_file:
        fasta_file.write(fasta_content)

    # remove the subset gfa files
    for ARO, path in congruent_paths:
        subset_gfa_file = gfa_input.split(".")[0] + "_" + ARO + "_subset.gfa"
        try:
            os.remove(subset_gfa_file)
        except OSError:
            print("Error: %s : %s" % (subset_gfa_file, "File not found"))

if __name__ == '__main__':
    graphaligner_output = sys.argv[1]
    gfa_input = sys.argv[3]
    summarized_hits = process_graphaligner_output(graphaligner_output)
    # must make ARO column unique. Add _instance_1,2,3 etc as a suffix
    summarized_hits["ARO"] = summarized_hits.groupby("ARO").cumcount().add(1).astype(str).radd(summarized_hits["ARO"] + "inst")
    summarized_hits.to_csv(sys.argv[2], sep="\t", index=False)
    get_path_sequences(summarized_hits, gfa_input)