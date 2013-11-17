#!/usr/bin/env python
import sys
from util import motif_finder

if __name__ == "__main__":
    f = open(sys.argv[1])
    lines = f.readlines()
    k = int(lines[0].strip())
    dna = motif_finder.MotifFinder(lines[1].strip())
    other_dnas = list()
    for i in range(2, len(lines)):
        other_dnas.append(motif_finder.MotifFinder(lines[i].strip()))
    median_motif = dna.median_kmer(other_dnas, k)
    print median_motif

