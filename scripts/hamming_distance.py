#!/usr/bin/env python
import sys

from util import dna_stats

if __name__ == "__main__":
    f = open(sys.argv[1])
    l = map(lambda x: x.strip(), f.readlines())
    d = dna_stats.DNAMultiStats(l)
    print d.hamming_distance()
