#!/usr/bin/env python
from util.euler_util import EulerUtil
import sys

if __name__ == "__main__":
    adj_list = dict()
    f = open(sys.argv[1])
    d = int(f.readline().strip())
    lines = map(lambda a: a.strip(), f.readlines())
    for l in lines:
        pair = l.split("|")
	k = len(pair[0])
        prefix = (pair[0][:-1], pair[1][:-1])
        suffix = (pair[0][1:], pair[1][1:])
        try:
            adj_list[prefix] += [ suffix ]
        except:
            adj_list[prefix] = [ suffix ]
    
    e = EulerUtil(adj_list)
    ret_list = e.euler_path()
    orig_first_ret_list = map(lambda a: a[0], ret_list)
    orig_second_ret_list = map(lambda a: a[1], ret_list)
    orig_first_string = e.reconstruct_string(orig_first_ret_list)
    orig_second_string = e.reconstruct_string(orig_second_ret_list)
    print orig_first_string + orig_second_string[-d-k:]
