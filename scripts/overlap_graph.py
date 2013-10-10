#!/usr/bin/env python
import sys

from util import dna_transformer
from util import parserUtil

if __name__ == "__main__":
    f = open(sys.argv[1])
    l = map(lambda a: a.strip(), f.readlines()) 
    parser = parserUtil.fastaParser(l)
    parser.parse()
    s = dna_transformer.DNAMultiTransformer(parser.sequences, parser.names)
    overlap_dict = s.k_overlap(3)
    for k in overlap_dict.keys():
        for v in overlap_dict[k]:
            print k, v
