#!/usr/bin/env python
import sys

from util import protein_util, protein_transformer

if __name__ == "__main__":
    f = open(sys.argv[1])
    p = protein_transformer.ProteinTransformer("")
    l = f.readline().strip()
    spectrum = map(lambda a: int(a), l.split())
    possible_peptides = p.cyclopeptide_sequence(spectrum)
    mass_dict = dict()
    for peptide in possible_peptides:
        mass = map(lambda a: str(protein_util.protein_mass_int[a]), peptide)
        mass_dict["-".join(mass)] = 0
    for k in sorted(mass_dict.keys(), reverse=True):
        print k,
        
