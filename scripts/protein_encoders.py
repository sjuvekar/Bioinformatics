#!/usr/bin/env python
import sys

from util import rna_transformer, dna_transformer

if __name__ == "__main__":
    f = open(sys.argv[1])
    d = f.readline().strip()
    dna = dna_transformer.DNATransformer(d)
    rna = rna_transformer.RNATransformer(dna.transcribe())
    p = f.readline().strip()
    encoders = rna.encoding_strings(p)
    for encoder in encoders:
        print rna_transformer.RNATransformer(encoder).reverse_transcribe(),
