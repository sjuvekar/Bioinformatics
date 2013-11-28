#!/usr/bin/env python
import sys
from util import motif_finder

if __name__ == "__main__":
    f = open(sys.argv[1])
    lines = f.readlines()
    nums = lines[0].strip().split()
    k = int(nums[0])
    t = int(nums[1])
    N = int(nums[1])
    dna = motif_finder.MotifFinder(lines[1].strip())
    other_dnas = list()
    for i in range(2, len(lines)):
        other_dnas.append(motif_finder.MotifFinder(lines[i].strip()))
    
    d = dict()
    for i in range(20):
	print i
	best_motif = dna.gibbs_sampler_search(other_dnas, k, N)
	try:
		d[" ".join(best_motif)] += 1
	except:
		d[" ".join(best_motif)] = 1

    best_d = max(d, key=d.get)
    print d[best_d]
    for x in best_d.split():
	print x

