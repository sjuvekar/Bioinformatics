#!/usr/bin/env python
import sys

from util import protein_transformer

if __name__ == "__main__":
    f = open(sys.argv[1])
    l = f.readline().strip()
    p = protein_transformer.ProteinTransformer(l)
    spectrum = p.cyclospectrum()
    for s in spectrum:
        print s,
