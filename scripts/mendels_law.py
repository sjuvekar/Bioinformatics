#!/usr/bin/env python
"""
K with XX alleles, M with Xy alleles and N with y alleales
"""
import sys
import pandas

def dominant_probability(k, m, n):
    """
    This is 3 X 3 matrix depicting the probability of dominant trait for (K, M, N) rows and columns
    """
    trait_probabilities = pandas.DataFrame([[1., 1., 1.], [1., 0.75, 0.5], [1., 0.5, 0.]])
    """ 
    This is selection probability
    """
    selection_probabilities = pandas.DataFrame([[float(k*(k-1)/2.), float(k*m), float(k*n)], 
                                                [float(m*k), float(m*(m-1)/2.), float(m*n)], 
                                                [float(n*k), float(n*m), float(n*(n-1)/2.)]])
    selection_probabilities = selection_probabilities / float( (k+m+n) * (k+m+n-1) / 2)
    prod = trait_probabilities * selection_probabilities
    return prod[0][0] + prod[1][0] + prod[1][1] + prod[2][0] + prod[2][1] + prod[2][2]

if __name__ == "__main__":
    print dominant_probability(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))
    
