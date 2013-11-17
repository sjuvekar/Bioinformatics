#!/usr/bin/env python
import sys
import pandas
from util import motif_finder

if __name__ == "__main__":
    f = open(sys.argv[1])
    lines = f.readlines()
    dna = motif_finder.MotifFinder(lines[0].strip())
    k = int(lines[1].strip())
    chars = lines[2].split()
    profile_list = list()
    for i in range(3, len(lines)):
        probs = map(lambda a: float(a), lines[i].strip().split())
        profile_list.append(probs)
    profile = pandas.DataFrame(profile_list, columns = chars)
    best_kmer = dna.profile_most_probable_kmer(k, profile)
    print best_kmer

