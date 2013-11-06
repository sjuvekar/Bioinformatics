#!/usr/bin/env python
import sys
from util import oriC_finder

if __name__ == "__main__":
    f = open(sys.argv[1])
    line = f.readline().strip().split()
    seq = line[0]
    kmer_len = int(line[1])
    thresh = int(line[2])
    oriC = oriC_finder.OriCFinder(seq)
    kmer_freq = oriC.frequent_words_with_mismatch(kmer_len, thresh)
    max_kmer_freq = max(kmer_freq.values())
    for k in kmer_freq.keys():
        if kmer_freq[k] == max_kmer_freq:
            print k,

