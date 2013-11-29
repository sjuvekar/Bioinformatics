#!/usr/bin/env python
from util.euler_util import EulerUtil
from util.string_mutator import StringMutator
import sys

if __name__ == "__main__":
    adj_list = dict()
    f = open(sys.argv[1])
    k = int(f.readline().strip())
    s = StringMutator()
    bases = ["0", "1"]
    for kmer in s.lengthKKmers(k-1, bases):
        suffix = kmer[1:]
        neighbors = map(lambda a: suffix + a, bases)
        adj_list[kmer] = neighbors
    
    e = EulerUtil(adj_list)
    ret_list = e.euler_tour()
    orig_string = e.reconstruct_string(ret_list)
    print orig_string[:1-k]
