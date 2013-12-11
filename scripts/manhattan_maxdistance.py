#!/usr/bin/env python
import sys
import pandas

def manhattan_maxdistance(rows, cols, down_matrix, right_matrix):
    optimal_distance = pandas.DataFrame(columns=range(cols+1), index=range(rows+1))
    optimal_distance[cols][rows] = 0

    for row in range(rows, -1, -1):
        for col in range(cols, -1, -1):
            down_distance = 0
            right_distance = 0
            if col < cols:
                right_distance = right_matrix[col][row] + optimal_distance[col+1][row]
            if row < rows:
                down_distance = down_matrix[col][row] + optimal_distance[col][row+1]
            optimal_distance[col][row] = max(right_distance, down_distance)

    return optimal_distance[0][0]


if __name__ == "__main__":
    f = open(sys.argv[1])
    rows = int(f.readline())
    cols = int(f.readline())
    down_matrix = pandas.DataFrame(columns=range(cols+1), index=range(rows))
    right_matrix = pandas.DataFrame(columns=range(cols), index=range(rows+1))
    for i in range(rows):
        next_line = f.readline().strip()
        down_matrix.ix[i] = map(lambda a:int(a), next_line.split())
    dash = f.readline().strip()

    assert dash == "-"

    for i in range(rows+1):
        next_line = f.readline().strip()
        right_matrix.ix[i] = map(lambda a:int(a), next_line.split())

    print manhattan_maxdistance(rows, cols, down_matrix, right_matrix)
    

