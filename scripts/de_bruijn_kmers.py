#!/usr/bin/env python
import sys

if __name__ == "__main__":
    f = open(sys.argv[1])
    lines = map(lambda a: a.strip(), f.readlines())
    kmer_dict = dict()
    for l in lines:
        suffix = l[1:]
        prefix = l[:-1]
        try:
            kmer_dict[prefix][suffix] = 0
        except:
            kmer_dict[prefix] = dict()
            kmer_dict[prefix][suffix] = 0

    for k in sorted(kmer_dict.keys()):
        print k, "->", ",".join(sorted(kmer_dict[k].keys()))


