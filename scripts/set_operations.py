#!/usr/bin/env python 
import sys

def parse(s):
    x = s.split(", ")
    x[0] = x[0][1:]
    x[-1] = x[-1][:-1]
    return set(x)

def pprint(s):
    print "{" + ", ".join(s) + "}"

if __name__ == "__main__":
    f = open(sys.argv[1])
    l = f.readlines()
    m = int(l[0].strip())
    complete_set = set(map(lambda a: str(a), range(1, m+1)))
    x = parse(l[1].strip())
    y = parse(l[2].strip())
    
    pprint(x.union(y))
    pprint(x.intersection(y))
    pprint(x.difference(y))
    pprint(y.difference(x))
    pprint(complete_set.difference(x))
    pprint(complete_set.difference(y))
