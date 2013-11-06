#!/usr/bin/env python
import sys
from util import oriC_finder

if __name__ == "__main__":
    f = open(sys.argv[1])
    seq = f.readline().strip()
    line = f.readline().strip().split()
    kmer_len = int(line[0])
    thresh = int(line[1])
    oriC = oriC_finder.OriCFinder(seq)
    kmer_freq = oriC.frequent_words_with_mismatch_and_reverse_complement(kmer_len, thresh)
    max_kmer_freq = max(kmer_freq.values())
    for k in kmer_freq.keys():
        if kmer_freq[k] == max_kmer_freq:
            print k,

