#!/usr/bin/env python
from util.euler_util import EulerUtil
import sys

if __name__ == "__main__":
    adj_list = dict()
    f = open(sys.argv[1])
    lines = map(lambda a: a.strip(), f.readlines())
    for l in lines:
        prefix = l[:-1]
        suffix = l[1:]
        try:
            adj_list[prefix] += [ suffix ]
        except:
            adj_list[prefix] = [ suffix ]
    
    e = EulerUtil(adj_list)
    ret_list = e.all_contigs()
    for r in ret_list:
	new_r = map(lambda a: a[0], r)
        print e.reconstruct_string(new_r)
