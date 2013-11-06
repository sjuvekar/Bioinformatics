#!/usr/bin/env python
import sys
from util import oriC_finder

if __name__ == "__main__":
    f = open(sys.argv[1])
    seq = f.readline().strip()
    oriC = oriC_finder.OriCFinder(seq)
    skew_positions = oriC.min_gc_skew()
    for p in skew_positions:
        print p+1,

