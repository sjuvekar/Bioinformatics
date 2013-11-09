#!/usr/bin/env python
import sys

from util import protein_transformer

if __name__ == "__main__":
    p = protein_transformer.ProteinTransformer("")
    f = open(sys.argv[1])
    l = f.readline().strip()
    spectrum = map(lambda a: int(a), l.split())
    spectrum = sorted(spectrum)
    convolution_spectrum = p.convolution_spectrum(spectrum)
    for c in sorted(convolution_spectrum):
        print c,
