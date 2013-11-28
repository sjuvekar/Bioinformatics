#!/usr/bin/env python
import sys
from util import string_mutator
from util import dna_graph_util

if __name__ == "__main__":
    f = open(sys.argv[1])
    lines = f.readlines()
    k = int(lines[0].strip())-1
    string = lines[1].strip()
    s = string_mutator.StringMutator()
    compositions = s.lexicographic_kmers(string, k)
    graph_util = dna_graph_util.DNAGraphUtil(list(set(compositions)))
    adj_list = graph_util.adjacency_list()
    for a in sorted(adj_list.keys()):
        print a, "->", ",".join(adj_list[a])

