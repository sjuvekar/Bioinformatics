#!/usr/bin/env python
import sys
from util import oriC_finder

if __name__ == "__main__":
    f = open(sys.argv[1])
    pattern = f.readline().strip()
    seq = f.readline().strip()
    oriC = oriC_finder.OriCFinder(seq)
    positions = oriC.pattern_positions(pattern)
    for p in positions:
        print p,
