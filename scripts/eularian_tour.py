#!/usr/bin/env python
import copy
import sys

adj_list = dict()

def approx_euler_starting_from(start_vertex):
    if len(adj_list) == 0:
        return list()
    bifurcation_vertex = None
    ret_list = [start_vertex]
    curr_vertex = start_vertex
    next_vertex = None
    while next_vertex != start_vertex:
        # Find neighbors
        neighbors = adj_list[curr_vertex]

        # Set up first neighbor as next
        next_vertex = neighbors[0]

        # Modify lists
        # 1) Add next_vertex to ret_list
        ret_list.append(next_vertex)
        # 2) Remove next vertex from neighbors. The edge should never be visited
        neighbors.remove(next_vertex)
        if len(neighbors) == 0:
            adj_list.pop(curr_vertex)

        # Set curr_vertex 
        curr_vertex = next_vertex

    return ret_list


def find_bifurcation_vertex(approx_tour):
    for r in approx_tour:
        if r in adj_list:
            return r
    return None


def merge_tours(initial_tour, new_tour, bifurcation_vertex):
    split_index = initial_tour.index(bifurcation_vertex)
    return initial_tour[:split_index+1] + new_tour[1:] + initial_tour[split_index+1:]


def approx_euler():
    starting_vertex = adj_list.keys()[0]
    ret_list = approx_euler_starting_from(starting_vertex)
    while len(adj_list.keys()) > 0:
        next_start = find_bifurcation_vertex(ret_list)
        if not next_start:
            print "No next start"
        recursed_list = approx_euler_starting_from(next_start)
        ret_list = merge_tours(ret_list, recursed_list, next_start)
    return ret_list
        

if __name__ == "__main__":
    f = open(sys.argv[1])
    lines = map(lambda a: a.strip(), f.readlines())
    for l in lines:
        tok = l.split(" -> ")
        k = int(tok[0])
        vals = map(lambda a: int(a), tok[1].split(","))
        adj_list[k] = vals
    ret_list = approx_euler()
    print "->".join(map(lambda a: str(a), ret_list))
