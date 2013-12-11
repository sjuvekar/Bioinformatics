from util.dna_transformer import DNAMultiTransformer
import numpy

class DNAGraphUtil(DNAMultiTransformer):

    def __init__(self, input_dnas, input_names=None):
        """
        Constructor keeps track of prefixes of the input dnas
        """
        DNAMultiTransformer.__init__(self, input_dnas, input_names)
        self.dna_prefixes = dict()
        for d in input_dnas:
            curr_prefix = d[:-1]
            try:
                self.dna_prefixes[curr_prefix][d] = 0
            except:
                self.dna_prefixes[curr_prefix] = dict()
		self.dna_prefixes[curr_prefix][d] = 0
            
        # For DAG longest path
        self.longest_path_neighbors = dict()

 
    def adjacency_list(self):
        adj_list = dict()
        for d in self.dna_transformers:
            dna = d.DNA
            dna_suffix = dna[1:]
	    if dna_suffix in self.dna_prefixes:
                for other_dna in self.dna_prefixes[dna_suffix]:
                    try:
                        adj_list[dna] += [other_dna]
                    except:
                        adj_list[dna] = [other_dna]			    
        return adj_list

    
    def dag_longest_path(self, adj_list, source, dest):
        """
        adj_list is a dict from node -> (node, distance)
        """
        if source == dest:
            self.longest_path_neighbors[dest] = (dest, 0)
            return
        if source not in adj_list.keys():
            self.longest_path_neighbors[source] = (source, -numpy.inf)
            return
        longest_path_distance = 0
        longest_path_neighbor = None
        for neighbor_distance in adj_list[source]:
            neighbor = neighbor_distance[0]
            distance = neighbor_distance[1]
            if neighbor not in self.longest_path_neighbors.keys():
                self.dag_longest_path(adj_list, neighbor, dest)
            longest_neighbor_distance = self.longest_path_neighbors[neighbor]
            best_distance = longest_neighbor_distance[1]
            if best_distance < 0:
                continue
            if best_distance + distance > longest_path_distance:
                longest_path_distance = best_distance + distance
                longest_path_neighbor = neighbor
                
        self.longest_path_neighbors[source] = (longest_path_neighbor, longest_path_distance)
