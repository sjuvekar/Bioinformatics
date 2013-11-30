import copy

class EulerUtil(object):

    def __init__(self, adj_list):
        """
        Constructor creates the adj_list which is modified in the class. It also saves a copy of original adj_list
        """
        self.adj_list = copy.deepcopy(adj_list)
        self.orig_adj_list = adj_list
        self.inverse_adj_list = dict()
        self.neighbor_offset = dict()

    def construct_inverse_adj_list(self):
        """
        Create inverse of adj_list by inverting edges
        """
        for a in self.adj_list.keys():
            for n in self.adj_list[a]:
                try:
                    self.inverse_adj_list[n] += [a]
                except:
                    self.inverse_adj_list[n] = [a]


    def approx_euler_starting_from(self, start_vertex, end_vertex = None):
        """
        Creates a subtour of Euler tour starting at start_vertex and ending at end_vertex. 
        If end_vertex is None, creates a sub-cycle ending at start_vertex
        """
        if len(self.adj_list) == 0:
            return list()
    
        if not end_vertex:
            end_vertex = start_vertex
        
        ret_list = [start_vertex]
        curr_vertex = start_vertex
        next_vertex = None
        while next_vertex != end_vertex and curr_vertex in self.adj_list.keys():
            # Find neighbors
            neighbors = self.adj_list[curr_vertex]

            # Set up first neighbor as next
            next_vertex = neighbors[0]

            # Modify lists
            # 1) Add next_vertex to ret_list
            ret_list.append(next_vertex)
            # 2) Remove next vertex from neighbors. The edge should never be visited
            neighbors.remove(next_vertex)
            if len(neighbors) == 0:
                self.adj_list.pop(curr_vertex)

            # Set curr_vertex 
            curr_vertex = next_vertex

        return ret_list


    def find_bifurcation_vertex(self, approx_tour):
        """
        Given an incomplete Euler tour, find a vertex that still has an unexplored neighbor
        """
        for r in approx_tour:
            if r in self.adj_list:
                return r
        return None


    def merge_tours(self, initial_tour, new_tour, bifurcation_vertex):
        """
        Merge a new Euler tour into old one
        """
        split_index = initial_tour.index(bifurcation_vertex)
        return initial_tour[:split_index+1] + new_tour[1:] + initial_tour[split_index+1:]


    def euler_tour(self):
        """
        Create an eulerian tour
        """
        starting_vertex = self.adj_list.keys()[0]
        ret_list = self.approx_euler_starting_from(starting_vertex)
        while len(self.adj_list.keys()) > 0:
            next_start = self.find_bifurcation_vertex(ret_list)
            if not next_start:
                print "No next start"
            recursed_list = self.approx_euler_starting_from(next_start)
            ret_list = self.merge_tours(ret_list, recursed_list, next_start)
        return ret_list
        

    def find_euler_start_end(self):
        """
        Find the start and end vertex of Euler path. Start vertex has one more outgoing edge ans end vertex has one more incoming edge
        it stores this information in neighbor_offset dict
        """
        self.neighbor_offset = dict()
        for n in self.adj_list.keys():
            self.neighbor_offset[n] = len(self.adj_list[n])
        for n in self.adj_list.keys():
            for a in self.adj_list[n]:
                if a not in self.neighbor_offset.keys():
                    self.neighbor_offset[a] = 0
                self.neighbor_offset[a] -= 1
        start_vertex = None
        end_vertex = None
        for n in self.neighbor_offset.keys():
            if self.neighbor_offset[n] == 1:
                start_vertex = n
            elif self.neighbor_offset[n] == -1:
                end_vertex = n
        return (start_vertex, end_vertex)


    def euler_path(self):
        """
        Creates an Euler path. Note that Euler path has a unique start and end vertex and in this case
        has to be found first
        """
        (start_vertex, end_vertex) = self.find_euler_start_end()
        base_path = self.approx_euler_starting_from(start_vertex, end_vertex)
        while end_vertex in self.adj_list.keys():
            new_path = self.approx_euler_starting_from(end_vertex)
            base_path += new_path[1:]
        while len(self.adj_list.keys()) > 0:
            next_start = self.find_bifurcation_vertex(base_path)
            if not next_start:
                print "No next start"
            recursed_list = self.approx_euler_starting_from(next_start)
            base_path = self.merge_tours(base_path, recursed_list, next_start)
        return base_path


    def reconstruct_string(self, path):
        """
        Reconstructs the strings from the path. assume that all labels on path are strings.
        """
        ans = path[0]
        for i in range(1, len(path)):
            ans += path[i][-1]
        return ans

    
    def build_contig_graph(self):
        """
        Builds a contig geaph with all vertices with more than one incoming and outgoing edges
        duplicated (in+ out) times
        """
        contig_graph = dict()
        self.construct_inverse_adj_list()
        vertex_multiplicities = dict()
        for a in self.adj_list.keys():
            fwd_neighbor_count = len(self.adj_list[a])
            vertex_multiplicities[a] = fwd_neighbor_count
        for i in self.inverse_adj_list.keys():
            inv_neighbor_count = len(self.inverse_adj_list[i])
            try:
                vertex_multiplicities[i] += inv_neighbor_count
            except:
                vertex_multiplicities[i] = inv_neighbor_count
            
        # Modify to remove simple vertices
        for v in vertex_multiplicities.keys():
            if (v not in self.adj_list.keys() or len(self.adj_list[v]) == 1) and (v not in self.inverse_adj_list.keys() or len(self.inverse_adj_list[v]) == 1):
                vertex_multiplicities[v] = 0

        # Build the graph now
        for a in self.adj_list.keys():
            for n in self.adj_list[a]:
                new_a = (a, vertex_multiplicities[a])
                new_n = (n, vertex_multiplicities[n])
                try:
                    contig_graph[new_a] += [ new_n ]
                except:
                    contig_graph[new_a] = [new_n]
                if vertex_multiplicities[n] > 0:
                    vertex_multiplicities[n] -= 1
                if vertex_multiplicities[a] > 0:
                    vertex_multiplicities[a] -= 1

        return contig_graph
    

    def all_contigs(self):
        """
        return all contigs for the graph
        """
        contigs = list()
        self.adj_list = self.build_contig_graph()
        while len(self.adj_list) > 0:
            (start_vertex, end_vertex) = self.find_euler_start_end()
            ret_list = self.approx_euler_starting_from(start_vertex, end_vertex)
            contigs.append(ret_list)
        return contigs
