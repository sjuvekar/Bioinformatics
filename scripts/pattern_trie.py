#!/usr/bin/env python
import sys
from util import pattern_util

if __name__ == "__main__":
    f = open(sys.argv[1])
    p = pattern_util.PrefixTrie()
    for l in f.readlines():
        p.add_sequence(l.strip())
    p.set_id(1)
    p.pprint()
