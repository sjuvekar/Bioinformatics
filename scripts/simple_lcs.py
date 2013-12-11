#!/usr/bin/env python
from util.lcs_util import LCSUtil
import sys

if __name__ == "__main__":
    f = open(sys.argv[1])
    dna = LCSUtil(f.readline().strip())
    other_dna = LCSUtil(f.readline().strip())
    print dna.lcs(other_dna)
