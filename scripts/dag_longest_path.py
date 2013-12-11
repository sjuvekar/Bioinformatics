#!/usr/bin/env python
import sys
import numpy
from util.dna_graph_util import DNAGraphUtil

def dag_longest_path(adj_list, source, dest):
    """
    adj_list is a dict from node -> (node, distance)
    """
    if source == dest:
        longest_path_neighbors[dest] = (dest, 0)
        return
    if source not in adj_list:
        longest_path_neighbors[source] = (source, -numpy.inf)
        return
    longest_path_distance = 0
    longest_path_neighbor = None
    for neighbor_distance in adj_list[source]:
        neighbor = neighbor_distance[0]
        distance = neighbor_distance[1]
        if neighbor not in longest_path_neighbors.keys():
            dag_longest_path(adj_list, neighbor, dest)
        longest_neighbor_distance = longest_path_neighbors[neighbor]
        best_neighbor = longest_neighbor_distance[0]
        best_distance = longest_neighbor_distance[1]
        if best_distance + distance > longest_path_distance:
            longest_path_distance = best_distance + distance
            longest_path_neighbor = best_neighbor

    longest_path_neighbors[source] = (longest_path_neighbor, longest_path_distance)


if __name__ == "__main__":
    f = open(sys.argv[1])
    lines = f.readlines()
    source = int(lines[0])
    dest = int(lines[1])
    
    adj_list = dict()
    
    for l in range(2, len(lines)):
        edge = lines[l].strip()
        edge_tok = edge.split("->")
        edge_src = int(edge_tok[0])
        edge_dest = map(lambda a:int(a), edge_tok[1].split(":"))
        try:
            adj_list[edge_src] += [edge_dest]
        except:
            adj_list[edge_src] = [edge_dest]

    graph_util = DNAGraphUtil([])
    graph_util.dag_longest_path(adj_list, source, dest)
    print graph_util.longest_path_neighbors[source][1]
