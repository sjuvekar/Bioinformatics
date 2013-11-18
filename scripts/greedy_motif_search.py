#!/usr/bin/env python
import sys
from util import motif_finder

if __name__ == "__main__":
    f = open(sys.argv[1])
    lines = f.readlines()
    nums = lines[0].strip().split()
    k = int(nums[0])
    t = int(nums[1])
    dna = motif_finder.MotifFinder(lines[1].strip())
    other_dnas = list()
    for i in range(2, len(lines)):
        other_dnas.append(motif_finder.MotifFinder(lines[i].strip()))
    best_motif = dna.greedy_motif_search(other_dnas, k, True)
    for b in best_motif:
        print b

