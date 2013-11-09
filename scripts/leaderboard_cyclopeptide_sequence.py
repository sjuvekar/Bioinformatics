#!/usr/bin/env python
import sys

from util import protein_util, protein_transformer

if __name__ == "__main__":
    f = open(sys.argv[1])
    p = protein_transformer.ProteinTransformer("")
    leaderboard_len = int(f.readline().strip())
    l = f.readline().strip()
    spectrum = map(lambda a: int(a), l.split())
    best_peptide = p.leaderboard_cyclopeptide_sequence(leaderboard_len, spectrum)
    mass = map(lambda a: str(protein_util.protein_mass_int[a]), best_peptide)
    print "-".join(mass)
