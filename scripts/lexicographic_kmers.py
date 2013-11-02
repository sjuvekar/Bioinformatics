#!/usr/bin/env python
import sys

def lex_kmers(chars, len):
    if len == 1:
        return chars
    kmers = [""] + lex_kmers(chars, len-1)
    ret = []
    for c in chars:
        for k in kmers:
            ret += [str(c) + str(k)]
    return ret

if __name__ == "__main__":
    f = open(sys.argv[1])
    chars = f.readline().strip().split()
    len = int(f.readline())
    ret = lex_kmers(chars, len)
    for c in ret:
        print c
    
