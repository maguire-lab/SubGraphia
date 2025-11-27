#!/usr/bin/env python3

# this is a python script to extract subgraphs from a graph. 
# It first loads a graph from a GFA file, 
# loads in nodes of interest e.g. graphaligner output, 
# then extracts the subgraph of interest.
## USAGE: python3 subgraph_extract.py <gfa_input> <graphaligner_output> <extraction_radius>

import sys
import subprocess
import pandas as pd

# define a function to parse a graphaligner output gaf file to get the nodes of interest list as a tuple with the gene name
def parse_ga_output(ga_output):
    ga_df = pd.read_csv(ga_output, sep="\t")
    paths = ga_df["path"]
    # remove first character from each string in the list
    paths = [node[1:] for node in paths]
    # remove > or < and replace with , in each string in the list
    paths = [node.replace(">", ",").replace("<", ",") for node in paths]
    #create tuple of ga_df["ARO"] and paths
    nodes_of_interest = list(zip(ga_df["ARO"], paths))
    return nodes_of_interest

# define a function to use gfatools as a subprocess to extract subgraph

def extract_subgraph(gfa_input, nodes_of_interest, extraction_radius):
    #loop though the nodes of interest tuple and extract the subgraph for each (could potentially be parallelized)
    for node in nodes_of_interest:
        subprocess.run(["gfatools", "view", "-l", node[1], "-r", extraction_radius, gfa_input], stdout=open(str(node[0]) + ".gfa", "w"))

if __name__ == "__main__":
    gfa_input = sys.argv[1]
    ga_output = sys.argv[2]
    extraction_radius = sys.argv[3]
    nodes_of_interest = parse_ga_output(ga_output)
    extract_subgraph(gfa_input, nodes_of_interest, extraction_radius)
