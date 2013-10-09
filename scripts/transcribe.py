#!/usr/bin/env python
import sys

from util import dna_transformer

if __name__ == "__main__":
    f = open(sys.argv[1])
    l = f.readline().strip()
    d = dna_transformer.DNATransformer(l)
    ret = d.transcribe()
    print ret
