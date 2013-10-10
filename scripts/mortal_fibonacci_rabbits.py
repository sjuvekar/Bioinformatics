#!/usr/bin/env python
"""
We maintain a list of all rabbit ages of size m. 0th entry being recently born rabbits
All rabbits from position 1 onward reproduce
"""
import sys

def recurrence_rabits(n, m):
    rabbit_numbers = []
    for i in range(m):
        rabbit_numbers.append(0)
    rabbit_numbers[0] = 1
    for i in range(1,n):
        new_rabbits = sum(rabbit_numbers[1:])
        rabbit_numbers = [new_rabbits] + rabbit_numbers[:-1]
    return sum(rabbit_numbers)

if __name__ == "__main__":
    print recurrence_rabits(int(sys.argv[1]), int(sys.argv[2]))
    
