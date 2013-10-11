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
    rna_seq = rna_transformer.RNATransformer(dna_seq.transcribe())
    reverse_dna_seq = dna_transformer.DNATransformer(dna_seq.reverse_complement(), dna_seq.name)
    reverse_rna_seq = rna_transformer.RNATransformer(reverse_dna_seq.transcribe())

    proteins = rna_seq.all_possible_translations()
    reverse_proteins = reverse_rna_seq.all_possible_translations()
    for p in set(proteins).union(set(reverse_proteins)):
        print p

