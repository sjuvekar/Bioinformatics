#!/usr/bin/env python
import sys
from util import oriC_finder

if __name__ == "__main__":
    f = open(sys.argv[1])
    seq = f.readline().strip()
    params = f.readline().strip()
    (kmer_len, window_len, freq) = map(lambda a:int(a), params.split())
    oriC = oriC_finder.OriCFinder(seq)
    freq_kmers_dict = oriC.find_clumps(kmer_len, window_len, freq)
    for k in freq_kmers_dict.keys():
        print k,

