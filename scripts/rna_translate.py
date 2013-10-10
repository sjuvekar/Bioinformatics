#!/usr/bin/env python
import sys

from util import rna_transformer

if __name__ == "__main__":
    f = open(sys.argv[1])
    l = f.readline().strip()
    d = rna_transformer.RNATransformer(l)
    print d.translate()
