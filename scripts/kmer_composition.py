#!/usr/bin/env python
import sys
from util import string_mutator

if __name__ == "__main__":
    f = open(sys.argv[1])
    lines = f.readlines()
    k = int(lines[0].strip())
    string = lines[1].strip()
    s = string_mutator.StringMutator()
    compositions = s.lexicographic_kmers(string, k)
    for c in compositions:
        print c

