#!/usr/bin/env python
"""
The recurrence relation is f(n) = f(n-1) + kf(n-2)
"""
import sys

def recurrence_rabits(n, k):
    if k <= 2:
        return 1
    f_n_2 = 1
    f_n_1 = 1
    for idx in range(3, n+1):
        temp = f_n_1 + k * f_n_2
        f_n_2 = f_n_1
        f_n_1 = temp
    return f_n_1

if __name__ == "__main__":
    print recurrence_rabits(int(sys.argv[1]), int(sys.argv[2]))
    
