#!/usr/bin/env python3

# script to convert a path through a graph, for example produced by path_walk.py, into a fasta file

import networkx as nx
import gfapy
import sys
import csv
from Bio.Seq import Seq

def parse_to_directed_graph(gfa_input):
  gfa_graph = gfapy.Gfa.from_file(gfa_input)
  directed_graph = nx.DiGraph()

  ### create nodes
  for segment in gfa_graph.segments:
    directed_graph.add_node(segment.name, name = segment.name, sequence = segment.sequence, KC = segment.KC)

  #get edges from gfa_graph
  edges = [line for line in gfa_graph.lines if isinstance(line, gfapy.Line) and line.record_type == "L"]

  for edge in edges:
    directed_graph.add_edge(edge.from_segment.name, edge.to_segment.name, from_orient = edge.from_orient, to_orient = edge.to_orient)
  
  print("Graph loaded with " + str(directed_graph.number_of_nodes()) + " nodes and " + str(directed_graph.number_of_edges()) + " edges")
  return directed_graph

def path_to_fasta(directed_graph, path_file, overlap):
    """
    Function to convert an input path to a fasta file using the directed graph.
    Directed graph overlap must be supplied, intended for use with graph that have an overlap of k and a minimum segment length of k+1.
    """
    # read the path file and parse into a list of lists using csv
    with open(path_file, 'r') as f:
        reader = csv.reader(f, delimiter=',')
        paths = [list(map(str.strip, row)) for row in reader]

    # output a fasta file for each path 
    for path in paths:
        # if path contains only one node, write the sequence of that node to a fasta file
        if len(path) == 1:
            node = path[0]
            node_data = directed_graph.nodes[node]
            fasta_header = ">" + gfa_file.split(".")[0] + "_path_" + str(paths.index(path) + 1)
            with open(fasta_header[1:] + ".fasta", 'w') as fasta_file:
                fasta_file.write(fasta_header + "\n")
                fasta_file.write(node_data['sequence'] + "\n")
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
                    trimmed_sequence = node_data['sequence'][:-overlap]
                    path_sequence = trimmed_sequence
                elif node_orient == "-": #reverse complement
                    RC_sequence = Seq(node_data['sequence']).reverse_complement()
                    trimmed_sequence = str(RC_sequence)[:-overlap]
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
                    trimmed_sequence = node_data['sequence'][:-overlap]
                    #if the path direction is "out" add trimmed sequence to the end of the path sequence
                    if path_dir == "out":
                        path_sequence += trimmed_sequence
                    #if the path direction is "in" add trimmed sequence to the start of the path sequence
                    elif path_dir == "in":
                        path_sequence = trimmed_sequence + path_sequence
                elif node_orient == "-":
                    RC_sequence = Seq(node_data['sequence']).reverse_complement()
                    trimmed_sequence = str(RC_sequence)[:-overlap]
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
                    trimmed_sequence = str(RC_sequence)[:-overlap]
                    #if the path direction is "out" add trimmed sequence to the end of the path sequence
                    if path_dir == "out":
                        path_sequence += trimmed_sequence
                    #if the path direction is "in" add trimmed sequence to the start of the path sequence
                    elif path_dir == "in":
                        path_sequence = trimmed_sequence + path_sequence
                elif node_orient == "-":
                    trimmed_sequence = node_data['sequence'][:-overlap]
                    #if the path direction is "out" add trimmed sequence to the end of the path sequence
                    if path_dir == "out":
                        path_sequence += trimmed_sequence
                    #if the path direction is "in" add trimmed sequence to the start of the path sequence
                    elif path_dir == "in":
                        path_sequence = trimmed_sequence + path_sequence

        # write the path sequence to a fasta file
        #create a fasta header based on input gfa basename and a path index
        fasta_header = ">" + gfa_file.split(".")[0] + "_path_" + str(paths.index(path) + 1)

        with open(fasta_header[1:] + ".fasta", 'w') as fasta_file:
            fasta_file.write(fasta_header + "\n")
            fasta_file.write(path_sequence + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python path2fasta.py <gfa_file> <path_file> <overlap>")
        sys.exit(1)

    gfa_file = sys.argv[1]
    path_file = sys.argv[2]
    overlap = int(sys.argv[3])

    directed_graph = parse_to_directed_graph(gfa_file)
    path_to_fasta(directed_graph, path_file, overlap)



