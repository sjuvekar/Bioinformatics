#!/usr/bin/env python
import sys

from util import dna_transformer
from util import rna_transformer
from util import parserUtil

if __name__ == "__main__":
    f = open(sys.argv[1])
    l = map(lambda a: a.strip(), f.readlines())
    parser = parserUtil.fastaParser(l)
    parser.parse()
    dna_seq = dna_transformer.DNATransformer(parser.sequences[0], parser.names[0])
    dna_spliced = dna_transformer.DNATransformer(dna_seq.splice(parser.sequences[1:]))
    rna_seq = rna_transformer.RNATransformer(dna_spliced.transcribe())
    print rna_seq.translate()
