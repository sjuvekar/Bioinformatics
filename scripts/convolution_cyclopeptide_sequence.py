#!/usr/bin/env python
import sys

from util import protein_util, protein_transformer

if __name__ == "__main__":
    f = open(sys.argv[1])
    p = protein_transformer.ProteinTransformer("")
    convolution_len = int(f.readline().strip())
    leaderboard_len = int(f.readline().strip())
    l = f.readline().strip()
    spectrum = map(lambda a: int(a), l.split())
    best_peptide = p.convolution_cyclopeptide_sequence(spectrum, convolution_len, leaderboard_len) 
    mass = map(lambda a: str(a), best_peptide)
    print "-".join(mass)
