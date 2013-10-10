#!/usr/bin/env python
import sys

from util import dna_util
from util import dna_stats
from util import parserUtil

if __name__ == "__main__":
    f = open(sys.argv[1])
    l = map(lambda a: a.strip(), f.readlines()) 
    parser = parserUtil.fastaParser(l)
    parser.parse()
    s = dna_stats.DNAMultiStats(parser.sequences, parser.names)
    (consensus_matrix, consensus_string) = s.consensus_dna()
    print consensus_string
    for k in dna_util.bases:
        print k+":", " ".join(map(lambda x: str(x), consensus_matrix[k]))

