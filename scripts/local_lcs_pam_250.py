#!/usr/bin/env python
from util.lcs_util import LCSUtil
import sys

if __name__ == "__main__":
    f = open(sys.argv[1])
    dna = LCSUtil(f.readline().strip())
    dna.parse_score_matrix("../matrices/PAM250.txt")
    other_dna = LCSUtil(f.readline().strip())
    (best_score, seq1, seq2) = dna.graph_based_local_alignment(other_dna, 5, True)
    print best_score
    print seq1
    print seq2
