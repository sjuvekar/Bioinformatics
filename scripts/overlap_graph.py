#!/usr/bin/env python
import sys
from util import dna_graph_util

if __name__ == "__main__":
    f = open(sys.argv[1])
    lines = map(lambda a: a.strip(), f.readlines())
    graph_util = dna_graph_util.DNAGraphUtil(sorted(lines))
    adj_list = graph_util.adjacency_list()
    for a in sorted(adj_list.keys()):
	for b in adj_list[a]:
            print a, "->", b

