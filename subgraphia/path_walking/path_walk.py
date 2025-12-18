#!/usr/bin/env python3

# A python script to talk paths through a subgraph and return the paths as a list of segments and a fasta file of the sequence. 

import networkx as nx
import gfapy
import sys
import csv
from ete4 import NCBITaxa
import pandas as pd
import statistics
from Bio.Seq import Seq

def parse_to_directed_graph(gfa_input):
  gfa_graph = gfapy.Gfa.from_file(gfa_input)
  directed_graph = nx.DiGraph()

  ### create nodes
  for segment in gfa_graph.segments:
    directed_graph.add_node(segment.name, name = segment.name, sequence = segment.sequence, KC = segment.KC, DP = segment.DP)

  #get edges from gfa_graph
  edges = [line for line in gfa_graph.lines if isinstance(line, gfapy.Line) and line.record_type == "L"]

  # get overlap length from edges
  overlap = int(str(edges[0].overlap)[:-1])

  for edge in edges:
    directed_graph.add_edge(edge.from_segment.name, edge.to_segment.name, from_orient = edge.from_orient, to_orient = edge.to_orient)
  
  print("Graph loaded with " + str(directed_graph.number_of_nodes()) + " nodes and " + str(directed_graph.number_of_edges()) + " edges")
  return directed_graph, overlap

def ensure_segment_traversal(path, anchor_node, directed_graph):
    """
    nx.all_simple_paths does not take into account that segments have ends corresponding with the 5' and 3' ends of the sequence.
    This function examines a simple path and ensures that each segments is being traversed by the edges rather than edges connecting to just one end of a segment.
    A simple path does not take into account edge flipping and other assembly graph properties. Simplifies the problem as paths should be flowing out or in. 
    
    :param path: A list of node IDs representing a path through the graph, e.g. output of nx.all_simple_paths()
    :param anchor_node: The node ID paths are being anchored to (starting or ending node)
    :param directed_graph: Parsed .GFA graph from parse_to_directed_graph()
    :return: Trimmed path that ensures segment traversal, whole path if all segments are traversed correctly
    """

    node_orients = {}
    # if path_orient == "out":
    # Iterate through the path
    for i in range(len(path) - 1):
        current_node = path[i]
        next_node = path[i + 1]
        edge = directed_graph.get_edge_data(current_node, next_node)
        # set node orients based on edge attributes
        from_node_orient = edge['from_orient']
        to_node_orient = edge['to_orient']
        # Assign orientations to node_orients, add orients rather than overwrite
        if current_node not in node_orients:
            node_orients[current_node] = set()
        if next_node not in node_orients:
            node_orients[next_node] = set()
        node_orients[current_node].add(from_node_orient)
        node_orients[next_node].add(to_node_orient)
    # if any nodes in node_orients have more than one orientation, trim the path at that node (inclusive)
    for i in range(len(path)):
        node = path[i]
        if len(node_orients[node]) > 1:
            # Trim the path at this node (inclusive), keep side with anchor_node
            if anchor_node == path[0]:
                trimmed_path = path[:i]
            elif anchor_node == path[-1]:
                trimmed_path = path[i+1:]
            else:
                raise ValueError("Anchor node must be either the start or end of the path.")
            return trimmed_path
    return path

def remove_nested_paths(paths):
    """
    Deterministically removes paths that are nested within longer paths.

    Strategy:
    - Remove exact duplicate paths (keep first occurrence).
    - Sort unique paths by descending length, tie-breaking lexicographically.
    - Keep a path if it is not a subset of any already-kept (longer or equal) path.
    - Return the kept paths in their original input order.

    Parameters:
    paths: list of paths (each path is a list of nodes)

    Returns:
    A list of paths with nested paths removed (deterministic).
    """
    if not paths:
        return []

    # Keep track of original indices and drop exact duplicates deterministically.
    indexed = [(i, list(p)) for i, p in enumerate(paths)]
    seen = set()
    unique_indexed = []
    for idx, p in indexed:
        tup = tuple(p)
        if tup in seen:
            continue
        seen.add(tup)
        unique_indexed.append((idx, p))

    # Sort by descending length, then lexicographically to have a deterministic selection order.
    sorted_indexed = sorted(unique_indexed, key=lambda ip: (-len(ip[1]), tuple(ip[1])))

    kept = []
    kept_sets = []
    for idx, p in sorted_indexed:
        pset = set(p)
        # If p is a subset of any already kept path, skip it.
        if any(pset.issubset(k) for k in kept_sets):
            continue
        kept.append((idx, p))
        kept_sets.append(pset)

    # Return kept paths in their original input order (deterministic).
    kept_sorted_by_idx = [p for _, p in sorted(kept, key=lambda x: x[0])]
    return kept_sorted_by_idx

def collect_orientations(directed_graph, edge_list, anchor_node):
    """
    Collects the orientations of the nodes in an edge list
    returns a dictionary of node IDs and orientations
    """
    # create a dictionary to store the orientations of the nodes
    orient_dict = {}
    # loop through the edges in the edge list
    for edge in edge_list:
        # if edges are in_edges, get the orientation of the to node
        if edge[1] == anchor_node:
            orient_dict[edge[0]] = directed_graph.edges[edge]["to_orient"]
        # if edges are out_edges, get the orientation of the from node
        elif edge[0] == anchor_node:
            orient_dict[edge[1]] = directed_graph.edges[edge]["from_orient"]
    return orient_dict

def get_amr_nodes(alignment_file):
    """
    Function to parse a graphaligner alignment file and return a list of node IDs that align to the references used. 
    """
    node_list = []
    # open alignment file with pandas
    alignments_df = pd.read_csv(alignment_file, sep="\t")
    # select path column
    path_column = alignments_df["path"]
    #split each path by ">" or "<" and add to node list
    for path in path_column:
        nodes = path.replace(">", " ").replace("<", " ").split()
        # if node not in list, add it
        for node in nodes:
            if node not in node_list:
                node_list.append(node)
    return node_list

def extract_paths(directed_graph, target_nodes):
    """
    Parmeters:
    directed_graph: output of parse_to_directed_graph, 
    parsed graph should be a subgraph of interest and not a whole assembly graph 
    target_nodes : path output from graphaligner in format >nodeID>nodeID etc. Path cannot be forked. 
    Returns:
    full_paths : A list of node IDs representing paths to and from the target nodes
    """
    ## parse target nodes into a tuple of edge direction and node ID, e.g. ('>', 'node1') or ('<', 'node2')
    # if target nodes is a list, convert target nodes list to a string, necessary for command line input, may not be if read from the graphaligner output file
    if isinstance(target_nodes, list):
        target_nodes = "".join(target_nodes)
    # substitute > for >, and < for <, then split by ,
    target_nodes = target_nodes.replace(">", ",>,")
    target_nodes = target_nodes.replace("<", ",<,")
    target_nodes = target_nodes.split(",")
    # remove empty strings
    target_nodes = [node for node in target_nodes if node != ""]

    # create a tuple from the target nodes list with the direction and node ID
    target_nodes = [(target_nodes[i], target_nodes[i+1]) for i in range(0, len(target_nodes), 2)]

    ### DM: ADD A CHECK TO SEE IF THE TARGET NODES FORM A LINEAR PATH, IF NOT, PRINT A MESSAGE AND EXIT THE FUNCTION ###

    # The longest target node will be set as the anchor node from which edges and paths will be found. Any paths that do not contain the complete target node path will be discarded.
    # find the longest target node
    max_length = 0
    for node in target_nodes:
        # get length of node sequence
        node_length = len(directed_graph.nodes[node[1]]["sequence"])
        # if the length of the node is greater than the max length, set the max length to the node length and set this node as the anchor node
        if node_length > max_length:
            max_length = node_length
            anchor_node = node[1]

    # examine the edges connecting in/out of the target nodes to determine the direction of the path
    # find in and out edges of the first target node
    anchor_node_in_edges = directed_graph.in_edges(anchor_node)
    anchor_node_out_edges = directed_graph.out_edges(anchor_node)

    ## Look for ancestors and decendants of the target node
    # create empty containers for ancestors and decendants
    ancestors = []
    decendants = []
    only_ins = False
    only_outs = False
    both_in_out = False

    # Look at in and out edges to see which direction we need to look 
    # if only in edges exist, look for ancestors
    if len(anchor_node_in_edges) > 0 and len(anchor_node_out_edges) == 0:
        # sort to make deterministic
        ancestors = sorted(nx.ancestors(directed_graph, anchor_node))
        only_ins = True
    # if only out edges exist, look for decendants
    elif len(anchor_node_in_edges) == 0 and len(anchor_node_out_edges) > 0:
        decendants = sorted(nx.descendants(directed_graph, anchor_node))
        only_outs = True   
    # if both in and out edges exist, look for ancestors and decendants
    elif len(anchor_node_in_edges) > 0 and len(anchor_node_out_edges) > 0:
        ancestors = sorted(nx.ancestors(directed_graph, anchor_node))
        decendants = sorted(nx.descendants(directed_graph, anchor_node))
        both_in_out = True
    # if no edges exist, print a message that the node is isolated
    else:
        print("Node " + anchor_node + " has no neighbours")

    ### find paths ###
    full_paths = []

    # if only in edges exist, look for paths through nodes in ancestors set
    if only_ins:
        in_paths = []
        # find paths between nodes in ancestors and the target node
        for ancestor in ancestors:
            in_paths.extend([path for path in nx.all_simple_paths(directed_graph, ancestor, anchor_node)])
        
        # confirm segment traversal for each path
        in_paths_validated = []
        for path in in_paths:
            validated_path = ensure_segment_traversal(path, anchor_node, directed_graph)
            in_paths_validated.append(validated_path)

        # remove paths that do not contain all the target nodes
        in_paths_validated = [path for path in in_paths_validated if set([node[1] for node in target_nodes]).issubset(set(path))]

        # remove paths that are nested within longer paths
        in_paths_validated = remove_nested_paths(in_paths_validated)

        # Need to join in_paths together IF the in_edges connecting them to the target node have opposite target node orientations
        # find the orientations of the in edges connecting the target node to the ancestors
        #init a dictionary to store the orientations of the target nodes in the in edges
        # sort edge list to make deterministic before collecting orientations
        target_orient_dict = collect_orientations(directed_graph, sorted(list(anchor_node_in_edges)), anchor_node)
        
        # If dictionary has a mix of + and - orientations, join the paths
        if len(set(target_orient_dict.values())) > 1:

            # bin paths into + and - paths using the target_orient_dict
            right_orient_paths = []
            left_orient_paths = []
            for path in in_paths_validated:
                for node in sorted(target_orient_dict):
                    if node in path:
                        if target_orient_dict[node] == "+":
                            right_orient_paths.append(path)
                        elif target_orient_dict[node] == "-":
                            left_orient_paths.append(path)
            # join right and left orient paths
            for right_path in right_orient_paths:
                #remove last (target node) node from right path
                right_path = right_path[:-1]
                for left_path in left_orient_paths:
                    # reverse left path
                    left_path = left_path[::-1]
                    # combine right and left paths
                    full_paths.append(right_path + left_path)
        # if all in edges have the same orientation, add paths to full paths
        else:
            full_paths.extend(in_paths_validated)

        #if no paths are found, add target nodes to full paths with a message that no paths were found
        if len(full_paths) == 0:
            # create a node list of target nodes
            target_node_list = [node[1] for node in target_nodes]
            full_paths.append(target_node_list)
                
    # if only out edges exist, look for paths through nodes in decendants set
    elif only_outs:
        out_paths = []
        # find paths between nodes in decendants and the target node
        for decendant in decendants:
            out_paths.extend([path for path in nx.all_simple_paths(directed_graph, anchor_node, decendant)])
        # confirm segment traversal for each path
        out_paths_validated = []
        for path in out_paths:
            validated_path = ensure_segment_traversal(path, anchor_node, directed_graph)
            out_paths_validated.append(validated_path)

        # remove paths that do not contain all the target nodes
        out_paths_validated = [path for path in out_paths_validated if set([node[1] for node in target_nodes]).issubset(set(path))]

        # remove paths that are nested within longer paths
        out_paths_validated = remove_nested_paths(out_paths_validated)
        
        # Need to join out_paths together IF the out_edges connecting them to the target node have opposite target node orientations
        # find the orientations of the out edges connecting the target node to the decendants
        #init a dictionary to store the orientations of the target nodes in the out edges
        target_orient_dict = collect_orientations(directed_graph, sorted(list(anchor_node_out_edges)), anchor_node)
        
        # If dictionary has a mix of + and - orientations, join the paths
        if len(set(target_orient_dict.values())) > 1:

            # bin paths into + and - paths using the target_orient_dict
            right_orient_paths = []
            left_orient_paths = []
            for path in out_paths_validated:
                for node in sorted(target_orient_dict):
                    if node in path:
                        if target_orient_dict[node] == "+":
                            right_orient_paths.append(path)
                        elif target_orient_dict[node] == "-":
                            left_orient_paths.append(path)
            # join right and left orient paths
            for left_path in left_orient_paths:
                #reverse left path
                left_path = left_path[::-1]
                for right_path in right_orient_paths:
                    #remove the first node (anchor node) from the right path
                    right_path = right_path[1:]
                    full_paths.append(left_path + right_path)

        # if all out edges have the same orientation, add paths to full paths
        else:
            full_paths.extend(out_paths_validated)

        #if no paths are found, add target nodes to full paths with a message that no paths were found
        if len(full_paths) == 0:
            # create a node list of target nodes
            target_node_list = [node[1] for node in target_nodes]
            full_paths.append(target_node_list)

    # if both in and out edges exist, look for paths through nodes in ancestors and decendants
    elif both_in_out:
        #first look for paths from ancestors to the target node
        in_paths = []
        # find paths between nodes in ancestors and the target node
        for ancestor in ancestors:
            in_paths.extend([path for path in nx.all_simple_paths(directed_graph, ancestor, anchor_node)])
        # confirm segment traversal for each path
        in_paths_validated = []
        for path in in_paths:
            validated_path = ensure_segment_traversal(path, anchor_node, directed_graph)
            in_paths_validated.append(validated_path)

        # remove nested paths
        in_paths_validated = remove_nested_paths(in_paths_validated)

        #collect in_path orientations (sort edge list for deterministic behaviour)
        in_orient_dict = collect_orientations(directed_graph, sorted(list(anchor_node_in_edges)), anchor_node)
        
        # bin paths into + and - paths using the in_orient_dict
        right_in_paths = []
        left_in_paths = []
        for path in in_paths_validated:
            for node in sorted(in_orient_dict):
                if node in path:
                    if in_orient_dict[node] == "+":
                        right_in_paths.append(path)
                    elif in_orient_dict[node] == "-":
                        left_in_paths.append(path)

        # look for paths from the target node to decendants
        out_paths = []
        # find paths between nodes in decendants and the target node
        for decendant in decendants:
            out_paths.extend([path for path in nx.all_simple_paths(directed_graph, anchor_node, decendant)])
        # confirm segment traversal for each path
        out_paths_validated = []
        for path in out_paths:
            validated_path = ensure_segment_traversal(path, anchor_node, directed_graph)
            out_paths_validated.append(validated_path)

        # remove nested paths
        out_paths_validated = remove_nested_paths(out_paths_validated)

        #collect out_path orientations (sort edge list for deterministic behaviour)
        out_orient_dict = collect_orientations(directed_graph, sorted(list(anchor_node_out_edges)), anchor_node)

        # bin paths into + and - paths using the out_orient_dict
        right_out_paths = []
        left_out_paths = []
        for path in out_paths_validated:
            for node in sorted(out_orient_dict):
                if node in path:
                    if out_orient_dict[node] == "+":
                        right_out_paths.append(path)
                    elif out_orient_dict[node] == "-":
                        left_out_paths.append(path)
        
        # now join the paths in the appropriate orientations:
        # all right in paths with all right out paths
        for right_in_path in right_in_paths:
            #remove the last node (anchor node) from the right in path
            right_in_path = right_in_path[:-1]
            for right_out_path in right_out_paths:
                full_paths.append(right_in_path + right_out_path)
        
        # all right in paths with all left in paths flipped
        for right_in_path in right_in_paths:
            #remove the last node (anchor node) from the right in path
            right_in_path = right_in_path[:-1]
            for left_in_path in left_in_paths:
                # reverse left in path
                left_in_path = left_in_path[::-1]
                full_paths.append(right_in_path + left_in_path)

        # all left in paths with all left out paths
        for left_in_path in left_in_paths:
            #remove the last node (anchor node) from the left in path
            left_in_path = left_in_path[:-1]
            for left_out_path in left_out_paths:
                full_paths.append(left_in_path + left_out_path)

        # all left out paths flipped with all right out paths
        for left_out_path in left_out_paths:
            #reverse left out path
            left_out_path = left_out_path[::-1]
            for right_out_path in right_out_paths:
                #remove the first node (anchor node) from the right out path
                right_out_path = right_out_path[1:]
                full_paths.append(left_out_path + right_out_path)

        #Remove paths that do not contain all the target nodes
        full_paths = [path for path in full_paths if set([node[1] for node in target_nodes]).issubset(set(path))]

        #if no paths are found, add target nodes to full paths with a message that no paths were found
        if len(full_paths) == 0:
            # create a node list of target nodes
            target_node_list = [node[1] for node in target_nodes]
            full_paths.append(target_node_list)
    
    return full_paths

def extend_paths(directed_graph, target_nodes, full_paths):
    """
    This function extends the paths produced by extract_paths. 
    This works by finding paths from the terminal nodes of the extracted paths and joining them with the extracted paths.
    Parameters:
    directed_graph: output of parse_to_directed_graph, parsed graph should be a subgraph of interest and not a whole assembly graph 
    target_nodes : path output from graphaligner in format >nodeID>nodeID etc. Path cannot be forked.
    full_paths : A list of node IDs representing paths to extend
    Returns:
    extended_paths : A list of node IDs representing paths to and from the target nodes extended as much as possible given the subgraph
    """

    extended_paths = []

    #get anchor node from target nodes
    ## parse target nodes into a tuple of edge direction and node ID, e.g. ('>', 'node1') or ('<', 'node2')
    # if target nodes is a list, convert target nodes list to a string, necessary for command line input, may not be if read from the graphaligner output file
    if isinstance(target_nodes, list):
        anchor_query = "".join(target_nodes)
    # substitute > for >, and < for <, then split by ,
    anchor_query = anchor_query.replace(">", ",>,")
    anchor_query = anchor_query.replace("<", ",<,")
    anchor_query = anchor_query.split(",")
    # remove empty strings
    anchor_query = [node for node in anchor_query if node != ""]

    # create a tuple from the target nodes list with the direction and node ID
    anchor_query = [(anchor_query[i], anchor_query[i+1]) for i in range(0, len(anchor_query), 2)]

    ### DM: ADD A CHECK TO SEE IF THE TARGET NODES FORM A LINEAR PATH, IF NOT, PRINT A MESSAGE AND EXIT THE FUNCTION ###


    # The longest target node will be set as the anchor node from which edges and paths will be found. Any paths that do not contain the complete target node path will be discarded.
    # find the longest target node
    max_length = 0
    for node in anchor_query:
        # get length of node sequence
        node_length = len(directed_graph.nodes[node[1]]["sequence"])
        # if the length of the node is greater than the max length, set the max length to the node length and set this node as the anchor node
        if node_length > max_length:
            max_length = node_length
            anchor_node = node[1]

    # step one: extend right
    right_extended_paths = []
    # create a list of nodes already searched
    searched_nodes = []

    for path in full_paths:
        path_extended = False
        last_node = path[-1]
        right_paths = []
        
        #if last node has not been searched
        if last_node not in searched_nodes:
            right_paths = extract_paths(directed_graph, [">" + last_node])
            searched_nodes.append(last_node)

            # remove right_paths that do not share more than one node with the path
            right_paths = [right_path for right_path in right_paths if len(set(right_path).intersection(set(path))) > 1]

            # connect the right paths to the path as appropriate
            for right_path in right_paths:
                #if nodes are repeating in the path, continue
                if len(right_path) != len(set(right_path)):
                    continue # nx doesn't support loops, OG path will be added to right_extended_paths if only looping extenstions exist

                #determine which side of the path contains new extended nodes and which side shares nodes with the path
                #find index of last node in right path
                right_path_index = right_path.index(last_node)
                #if the path to the right of the last_node shares more than one node with path AND does not contain any nodes absent from the path:
                if len(set(path) - (set(path) - set(right_path[right_path_index:]))) > 1 and set(right_path[right_path_index:]).issubset(set(path)):
                    #check to see if the left side does not share more than one node with path
                    if len(set(path) - (set(path) - set(right_path[:right_path_index + 1]))) == 1: 
                        #join left side of the right_path with path
                        new_nodes = right_path[:right_path_index] #doesn't include last node
                        #reverse order of new nodes
                        new_nodes = new_nodes[::-1]
                        #join new nodes with all paths that share the last node
                        for p in full_paths:
                            if last_node == p[-1]:
                                right_extended_paths.append(p + new_nodes)
                        # remove duplicates deterministically
                        right_extended_paths = [list(x) for x in sorted({tuple(i) for i in right_extended_paths})]
                        path_extended = True
                                
                    else:
                        continue # both sides have membership with path, skip it
                #if the path to the left of the last_node shares more than one node with path AND does not contain any nodes absent from the path:
                elif len(set(path) - (set(path) - set(right_path[:right_path_index + 1]))) > 1 and set(right_path[:right_path_index + 1]).issubset(set(path)):
                    #check to see if the right side does not share more than one node with path
                    if len(set(path) - (set(path) - set(right_path[right_path_index:]))) == 1:
                        #join right side of the right_path with path
                        new_nodes = right_path[right_path_index + 1:] #doesn't include last node
                        #join new nodes with all paths that share the last node
                        for p in full_paths:
                            if last_node == p[-1]:
                                right_extended_paths.append(p + new_nodes)
                        # remove duplicates deterministically
                        right_extended_paths = [list(x) for x in sorted({tuple(i) for i in right_extended_paths})]
                        path_extended = True
                    else:
                        continue # both sides have membership with path, skip it
           
            # if a path has not been extended, add it to the right_extended_paths
            if path_extended == False:
                right_extended_paths.append(path)
             
    # remove nested paths from right_extended_paths
    right_extended_paths = remove_nested_paths(right_extended_paths)

    # remove nodes in each right_extended_path to the left of the anchor node
    trimmed_right_extended_paths = []
    for right_extended_path in right_extended_paths:
        trimmed_right_extended_paths.append(right_extended_path[right_extended_path.index(anchor_node):])

    # deduplicate and sort trimmed_right_extended_paths deterministically
    trimmed_right_extended_paths = [list(x) for x in sorted({tuple(i) for i in trimmed_right_extended_paths})]

    # now extend to the left
    left_extended_paths = []
    # create a list of nodes already searched
    searched_nodes = []

    for path in full_paths:
        path_extended = False
        first_node = path[0]
        left_paths = []
        # if first_node == anchor_node:
        #     continue
        if first_node not in searched_nodes:
            left_paths = extract_paths(directed_graph, ["<" + first_node])
            searched_nodes.append(first_node)

            # remove left_paths that do not share more than one node with the path
            left_paths = [left_path for left_path in left_paths if len(set(left_path).intersection(set(path))) > 1]

            # connect the left paths to the path as appropriate
            for left_path in left_paths:
                #if nodes are repeating in the path, continue
                if len(left_path) != len(set(left_path)):
                    continue # nx doesn't support loops, OG path will be added to left_extended_paths if only looping extenstions exist

                #determine which side of the path contains new extended nodes and which side shares nodes with the path
                #find index of first node in left path
                left_path_index = left_path.index(first_node)
                #if the path to the right of the first_node shares more than one node with path AND does not contain any nodes absent from the path:
                if len(set(path) - (set(path) - set(left_path[left_path_index:]))) > 1 and set(left_path[left_path_index:]).issubset(set(path)):
                    #check to see if the left side does share more than one node with path
                    if len(set(path) - (set(path) - set(left_path[:left_path_index + 1]))) == 1:
                        #join left side of the left_path with path
                        new_nodes = left_path[:left_path_index] #doesn't include first node
                        #join new nodes with all paths that share the first node
                        for p in full_paths:
                            if first_node == p[0]:
                                left_extended_paths.append(new_nodes + p)
                        # remove duplicates deterministically
                        left_extended_paths = [list(x) for x in sorted({tuple(i) for i in left_extended_paths})]
                        path_extended = True
                    else:
                        continue # both sides have membership with path, skip it
                #if the path to the left of the first_node shares more than one node with path AND does not contain any nodes absent from the path:
                elif len(set(path) - (set(path) - set(left_path[:left_path_index + 1]))) > 1 and set(left_path[:left_path_index + 1]).issubset(set(path)):
                    #check to see if the right side shares more than one node with path
                    if len(set(path) - (set(path) - set(left_path[left_path_index:]))) == 1:
                        #join right side of the left_path with path
                        new_nodes = left_path[left_path_index + 1:] #doesn't include first node
                        #reverse order of new nodes
                        new_nodes = new_nodes[::-1]
                        #join new nodes with all paths that share the first node
                        for p in full_paths:
                            if first_node == p[0]:
                                left_extended_paths.append(new_nodes + p)
                        # remove duplicates deterministically
                        left_extended_paths = [list(x) for x in sorted({tuple(i) for i in left_extended_paths})]
                        path_extended = True 
                    else:
                        continue # both sides have membership with path, skip it

            # if an path has not been extended add it to the left_extended_paths
            if path_extended == False:
                left_extended_paths.append(path)

    # remove nested paths from left_extended_paths
    left_extended_paths = remove_nested_paths(left_extended_paths)    
        
    # remove nodes in each left_extended_path to the right of the anchor node
    trimmed_left_extended_paths = []
    for left_extended_path in left_extended_paths:
        trimmed_left_extended_paths.append(left_extended_path[:left_extended_path.index(anchor_node)])

    # deduplicate and sort trimmed_left_extended_paths deterministically
    trimmed_left_extended_paths = [list(x) for x in sorted({tuple(i) for i in trimmed_left_extended_paths})]

    # join left and right extended paths in all possible combinations (iterate in deterministic order)
    for left_path in trimmed_left_extended_paths:
        for right_path in trimmed_right_extended_paths:
            extended_paths.append(left_path + right_path)

    # remove nested paths from extended_paths
    extended_paths = remove_nested_paths(extended_paths)

    ### NOT NECESSARILY A GOOD THING TO DO ###
    # remove paths that repeats nodes
    extended_paths = [path for path in extended_paths if len(path) == len(set(path))]

    return extended_paths

def taxid_path_filter(kraken_out, full_paths, target_nodes, alignment_file):
    # define ncbi taxonomy database
    ncbi = NCBITaxa()

    # get_amr_nodes
    amr_nodes = get_amr_nodes(alignment_file)

    # parse kraken output tsv file
    kraken_df = pd.read_csv(kraken_out, sep="\t", header=None, names=["classification_status","nodeID", "taxID", "length","LCA"])  
    kraken_df = kraken_df.dropna()
    
    # create a dictionary of nodeIDs and taxIDs
    node_taxid_dict = dict(zip(kraken_df["nodeID"], kraken_df["taxID"]))
    # convert node_taxid_dict nodeIDs to strings
    node_taxid_dict = {str(k):v for k,v in node_taxid_dict.items()}

    #create a dictionary of nodeIDs and length 
    node_length_dict = dict(zip(kraken_df["nodeID"], kraken_df["length"]))
    # convert node_length_dict nodeIDs to strings and lengths to integers
    node_length_dict = {str(k):int(v) for k,v in node_length_dict.items()}

    #parse target nodes
    if isinstance(target_nodes, list):
        target_nodes = "".join(target_nodes)
    # substitute > for >, and < for <, then split by ,
    target_nodes = target_nodes.replace(">", ",>,")
    target_nodes = target_nodes.replace("<", ",<,")
    target_nodes = target_nodes.split(",")
    # remove empty strings
    target_nodes = [node for node in target_nodes if node != ""]

    # create a tuple from the target nodes list with the direction and node ID
    target_nodes = [(target_nodes[i], target_nodes[i+1]) for i in range(0, len(target_nodes), 2)]

    #define congruent paths
    congruent_paths = []

    # Loop though each path in full_paths
    for path in full_paths:
        # define path lineage
        path_lineage = []
        #loop though each segment in the path
        for segment in path:
            #if segment in target nodes and is less than 5000bp long, continue to the next segment (prevents very long target nodes from not getting speciated when they should be)
            if segment in [node[1] for node in target_nodes] and node_length_dict[segment] < 5000:
                continue
            # if segment is an AMR node and is less than 5000bp long, continue to the next segment
            if segment in amr_nodes and node_length_dict[segment] < 5000:
                continue
            #if the segment is shorter than 165 bp (3*k should be a variable in future), continue to the next segment
            if node_length_dict[segment] < 165:
                continue

            # get lineage of the segment using the node_taxid_dict
            node_lineage = ncbi.get_lineage(node_taxid_dict[segment])

            # if the path lineage is empty, set the node lineage as the path lineage and continue to the next segment
            if len(path_lineage) == 0:
                if node_lineage is not None:
                    path_lineage = node_lineage

            # if node_lineage is nonetype, continue to the next segment. Lineage can be empty. 
            if node_lineage is None:
                continue

            # If the shorter of the two lineages is a subset of the longer lineage, set the longer lineage as the path lineage and continue to the next segment
            if len(node_lineage) < len(path_lineage):
                if set(node_lineage).issubset(set(path_lineage)):
                    continue
                # Allows some flexibility, only when lineage lenth difference is 1
                elif set(node_lineage[:-1]).issubset(set(path_lineage)) and abs(len(node_lineage) - len(path_lineage)) < 2:
                    continue
                else:
                    break

            elif len(node_lineage) > len(path_lineage):
                if set(path_lineage).issubset(set(node_lineage)):
                    path_lineage = node_lineage
                    continue
                # Allows some flexibility 
                elif set(path_lineage[:-1]).issubset(set(node_lineage)) and abs(len(node_lineage) - len(path_lineage)) < 2:
                    continue
                else:
                    break

            # If the two lineages are the same length, if they are the same, continue to the next segment
            elif len(node_lineage) == len(path_lineage):
                if node_lineage == path_lineage:
                    continue
                # Allows for the case where nodes are getting called differently at the final taxonomic rank
                elif node_lineage[:-1] == path_lineage[:-1]:
                    continue
                # if none of the above conditions are met, path is not congruent, break the loop
                else:
                    break
        # If the loop did not break, the path is congruent, add it to the congruent paths list
        else:
            congruent_paths.append(path)

    #if no congruent paths are found, add target nodes to congruent paths with a message that no paths were found
    if len(congruent_paths) == 0:
        # create a node list of target nodes
        target_node_list = [node[1] for node in target_nodes]
        congruent_paths.append(target_node_list)
        print("No congruent paths found, target node(s) added to congruent paths")

    print("Number of taxonomically congruent paths found: " + str(len(congruent_paths)))        
    return congruent_paths    

def gather_metadata(directed_graph, kraken_out, congruent_paths):
    """
    This function gathers metadata for the congruent paths, useful in downstream annalyses
    inputs:
    directed_graph: output of parse_to_directed_graph, parsed graph should be a subgraph of interest and not a whole assembly graph
    kraken_out: path to the kraken output file
    congruent_paths: list of paths that are taxonomically congruent
    returns:
    A tsv file with the following columns:
    path_id, path_length, tax_id_lineage, median_coverage, stdev_coverage, length_bp
    """

    # define ncbi taxonomy database
    ncbi = NCBITaxa()

    # parse kraken output tsv file
    kraken_df = pd.read_csv(kraken_out, sep="\t", header=None, names=["classification_status","nodeID", "taxID", "length","LCA"])  
    kraken_df = kraken_df.dropna()
    # create a dictionary of nodeIDs and taxIDs
    node_taxid_dict = dict(zip(kraken_df["nodeID"], kraken_df["taxID"]))
    # convert node_taxid_dict nodeIDs to strings
    node_taxid_dict = {str(k):v for k,v in node_taxid_dict.items()} 
    # create a dictionary of nodeIDs and length
    node_length_dict = dict(zip(kraken_df["nodeID"], kraken_df["length"]))
    # convert node_length_dict nodeIDs to strings and lengths to integers
    node_length_dict = {str(k):int(v) for k,v in node_length_dict.items()}

    # #create dictionary of nodeIDs and coverage from the directed graph
    # node_coverage_dict = {}
    # for node in directed_graph.nodes(data=True):
    #     node_coverage_dict[node[0]] = node[1]["DP"]

    # # convert node_coverage_dict nodeIDs to strings and coverage to integers
    # node_coverage_dict = {str(k):int(v) for k,v in node_coverage_dict.items()}

    # create a list to store the metadata for each path
    metadata = []

    # loop through each path in congruent_paths
    for path in congruent_paths:
        # assign a path ID which is: gfa_file.split(".")[0] + "_path_" + str(paths.index(path) + 1)
        path_id = gfa_input.split(".")[0] + "_path_" + str(congruent_paths.index(path) + 1)
        # get the path length
        path_length = sum([node_length_dict[node] for node in path]) - (55 * (len(path)))
        # # get the median coverage of the path
        # path_coverage = [node_coverage_dict[node] for node in path]

        # # if path has just one node, set median and stdev coverage to NA
        # if len(path) == 1:
        #     median_coverage = "NA"
        #     stdev_coverage = "NA"
        # else:
        #     # calculate median and standard deviation of the coverage
        #     median_coverage = statistics.median(path_coverage)
        #     stdev_coverage = statistics.stdev(path_coverage)
        
        # get the longest taxID lineage for the path
        path_lineage = []
        for segment in path:
            # get lineage of the segment using the node_taxid_dict
            node_lineage = ncbi.get_lineage(node_taxid_dict[segment])
            # if node_lineage is NoneType, continue
            if node_lineage is None:
                continue
            # if node_lineage is longer than the current path lineage, set the path lineage to the node lineage
            if len(node_lineage) > len(path_lineage):
                path_lineage = node_lineage
            else: # path lineage is longer or equal to node lineage, continue
                continue
        #cut path_lineage off at the species rank if the lineages goes beyond
        ranks = ncbi.get_rank(path_lineage)
        if "species" in ranks.values():
            #get the taxid for species rank
            for taxid, rank in ranks.items():
                if rank == "species":
                    species_taxid = taxid
            #get index of species taxid in path_lineage
            species_index = path_lineage.index(species_taxid)
            #cut path_lineage off at species index + 1
            path_lineage = path_lineage[:species_index + 1]

        # convert path lineage to a string
        path_lineage = ",".join([str(taxid) for taxid in path_lineage])

        # get the taxa names for the path lineage
        path_lineage_names = []
        for taxid in path_lineage.split(","):
            path_lineage_names.append(ncbi.get_taxid_translator([int(taxid)])[int(taxid)])
        # convert path lineage names to a string
        path_lineage_names = ",".join(path_lineage_names)

        #get the rank of the last taxid
        LCA_rank = []
        path_lineage_int = [int(taxid) for taxid in path_lineage.split(",")]
        ranks = ncbi.get_rank(path_lineage_int)
        print(ranks)
        for taxid in path_lineage_int:
            LCA_rank.append(ranks[taxid])
        LCA_rank = LCA_rank[-1]
        
        # add the metadata for the path to the metadata list
        metadata.append([path_id, path_length, path_lineage, path_lineage_names, LCA_rank, 
                        #  median_coverage, stdev_coverage
                        ])
    # create a dataframe from the metadata list
    metadata_df = pd.DataFrame(metadata, columns=["path_id", "path_length", "tax_id_lineage", "tax_ID_names", "LCA_rank", 
                                                #   "median_coverage", "stdev_coverage"
                                                  ])
    return metadata_df

def path_to_fasta(directed_graph, congruent_paths, overlap):
    """
    Function to convert an input path to a fasta file using the directed graph.
    Directed graph overlap must be supplied, intended for use with graph that have an overlap of k and a minimum segment length of k+1.
    """
    #create an empty multifasta
    fasta_content = ""

    # output a fasta file for each path 
    for path in congruent_paths:
        # if path contains only one node, add the sequence of that node to fasta content
        if len(path) == 1:
            node = path[0]
            node_data = directed_graph.nodes[node]
            fasta_header = ">" + gfa_input.split(".")[0] + "_path_" + str(congruent_paths.index(path) + 1)
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

        #create a fasta header based on input gfa basename and a path index
        fasta_header = ">" + gfa_input.split(".")[0] + "_path_" + str(congruent_paths.index(path) + 1)
        #write the fasta header and sequence to the fasta content
        fasta_content += fasta_header + "\n" + path_sequence + "\n"
    
    # write the fasta content to a file
    with open(gfa_input.split(".")[0] + "_paths.fasta", 'w') as fasta_file:
        fasta_file.write(fasta_content)

if __name__ == '__main__':
    gfa_input = sys.argv[1]
    target_nodes = sys.argv[2].split(",")
    kraken_out = sys.argv[3]
    alignment_file = sys.argv[4]
    
    directed_graph, overlap = parse_to_directed_graph(gfa_input)
    
    full_paths = extract_paths(directed_graph, target_nodes)
    print("Number of original paths found: " + str(len(full_paths)))
   
    # #run extend_paths until no new paths are found and the length of the paths does not change
    previous_paths = []
    extended_paths = full_paths
    
    count = 0

    while count < 10 and len(extended_paths) < 300:
        previous_paths = extended_paths
        count += 1
        print("Iteration: " + str(count))
        print("extending " + str(len(previous_paths)) + " paths")
        extended_paths = extend_paths(directed_graph, target_nodes, previous_paths)

    else:
        if len(extended_paths) > 0:
            print("Extending paths complete")
        else:
            #set full paths to extended paths
            extended_paths = full_paths
            print("Paths could not be extended, using original paths")

    print("Number of extended paths found: " + str(len(extended_paths)))

    congruent_paths = taxid_path_filter(kraken_out, extended_paths, target_nodes, alignment_file)
    with open(gfa_input.split(".")[0] + "_congruent_paths.csv", "w") as f:
        writer = csv.writer(f)
        writer.writerows(congruent_paths)

    metadata_df = gather_metadata(directed_graph, kraken_out, congruent_paths)
    metadata_df.to_csv(gfa_input.split(".")[0] + "_metadata.tsv", sep="\t", index=False)

    path_to_fasta(directed_graph, congruent_paths, overlap)
    
