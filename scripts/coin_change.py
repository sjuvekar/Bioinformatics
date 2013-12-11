#!/usr/bin/env python
import sys
import numpy

def coin_change(amount, changes):
    optimal_changes = dict()
    optimal_changes[0] = 0
 
    for money in range(1, amount+1):
        optimal_changes[money] = numpy.inf
        
        for change in changes:
            diff = money - change
            if diff >= 0:
                curr_optimal = optimal_changes[diff] + 1
                if curr_optimal < optimal_changes[money]:
                    optimal_changes[money] = curr_optimal
    return optimal_changes[amount]

if __name__ == "__main__":
    f = open(sys.argv[1])
    amount = int(f.readline())
    changes = map(lambda a: int(a), f.readline().split(","))
    print coin_change(amount, changes)
