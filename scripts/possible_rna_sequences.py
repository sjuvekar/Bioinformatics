#!/usr/bin/env python
import sys

from util import protein_util

if __name__ == "__main__":
    f = open(sys.argv[1])
    l = f.readline().strip()
    p = protein_util.ProteinUtil(l)
    print p.possible_rna_sequences()
