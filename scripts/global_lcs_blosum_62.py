#!/usr/bin/env python
from util.lcs_util import LCSUtil
import sys

if __name__ == "__main__":
    f = open(sys.argv[1])
    dna = LCSUtil(f.readline().strip())
    dna.parse_score_matrix("../matrices/Blosum62.txt")
    other_dna = LCSUtil(f.readline().strip())
    seq = dna.lcs(other_dna, 5, True)
    (seq1, seq2) = dna.construct_insdel_sequences(other_dna)
    print dna.matching_matrix[len(dna.DNA)][len(other_dna.DNA)]
    print seq1
    print seq2
