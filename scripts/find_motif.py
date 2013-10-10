#!/usr/bin/env python
import sys

from util import dna_stats

if __name__ == "__main__":
    f = open(sys.argv[1])
    l = f.readline().strip()
    text = f.readline().strip()
    d = dna_stats.DNAStats(l)
    ret = d.find_motif(text)
    print " ".join(map(lambda a:str(a), ret))
