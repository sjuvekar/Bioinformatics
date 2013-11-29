#!/usr/bin/env python
from util.euler_util import EulerUtil
import sys

if __name__ == "__main__":
    adj_list = dict()
    f = open(sys.argv[1])
    lines = map(lambda a: a.strip(), f.readlines())
    for l in lines:
        tok = l.split(" -> ")
        k = int(tok[0])
        vals = map(lambda a: int(a), tok[1].split(","))
        adj_list[k] = vals
    
    e = EulerUtil(adj_list)
    ret_list = e.euler_path()
    print "->".join(map(lambda a: str(a), ret_list))
