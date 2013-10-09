#!/usr/bin/env python
import sys

from util import dna_stats

if __name__ == "__main__":
    f = open(sys.argv[1])
    l = f.readline().strip()
    d = dna_stats.DNAStats(l)
    ret = d.count_bases()
    print ret["A"], ret["C"], ret["G"], ret["T"]
