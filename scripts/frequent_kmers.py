#!/usr/bin/env python
import sys
from util import oriC_finder

if __name__ == "__main__":
    f = open(sys.argv[1])
    seq = f.readline().strip()
    k = int(f.readline())
    oriC = oriC_finder.OriCFinder(seq)
    kmers = oriC.frequent_kmers(k)
    for kmer in kmers:
        print kmer,
    
