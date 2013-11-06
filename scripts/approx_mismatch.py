#!/usr/bin/env python
import sys
from util import oriC_finder

if __name__ == "__main__":
    f = open(sys.argv[1])
    pattern = f.readline().strip()
    seq = f.readline().strip()
    thresh = int(f.readline())
    oriC = oriC_finder.OriCFinder(seq)
    positions = oriC.approx_matches(pattern, thresh)
    for p in positions:
        print p,

