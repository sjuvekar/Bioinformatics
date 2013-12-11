from dna_transformer import DNATransformer
from dna_graph_util import DNAGraphUtil
import pandas

class LCSUtil(DNATransformer):
    """
    Class to compute the longest common subsequence of two 
    sequences
    """

    def __init__(self, input_dna, fasta_name = "Rosanlind_"):
        DNATransformer.__init__(self, input_dna, fasta_name)
        self.matching_matrix = pandas.DataFrame()
        self.backtrack_matrix = pandas.DataFrame()
        self.score_matrix = pandas.DataFrame()
        self.TOP = 0
        self.LEFT = 1
        self.DIAGONAL = 2
        
    
    def parse_score_matrix(self, score_matrix_path):
        """
        Score matrix from file
        """
        f = open(score_matrix_path)
        lines = f.readlines()
        bases = lines[0].strip().split()
        self.score_matrix = pandas.DataFrame(columns=bases, index=bases)
        for l in range(1, len(lines)):
            tok = lines[l].strip().split()
            ind = tok[0]
            values = map(lambda a:int(a), tok[1:])
            self.score_matrix.ix[ind] = values
        f.close()


    def lcs(self, other_dna, penalty = 0, use_score_matrix = False):
        """
        Finds lcs of two DNAs. first initializes matching and backtracking matrics to empty and recursively fills them
        """
        self.matching_matrix = pandas.DataFrame(columns=range(len(self.DNA)+1), index=range(len(other_dna.DNA)+1))
        self.backtrack_matrix = pandas.DataFrame(columns=range(len(self.DNA)+1), index=range(len(other_dna.DNA)+1))

        # Init
        for i in range(len(self.DNA)+1):
            self.matching_matrix[i][0] = -penalty * i
        for i in range(len(other_dna.DNA)+1):
            self.matching_matrix[0][i] = -penalty * i
        
        # Dynamic progamming
        self.lcs_rec_fill_matrices(other_dna, penalty, use_score_matrix)
        return self.construct_subsequence(other_dna)


    def lcs_rec_fill_matrices(self, other_dna, penalty = 0, use_score_matrix = False):
        """
        Fills the backtrack and matching matrices using Dynamic programming
        """
        for i in range(1, len(self.DNA)+1):
            for j in range(1, len(other_dna.DNA)+1):
                diag_best = 0
                left_best = 0
                top_best = 0
                self_char = self.DNA[i-1]
                other_char = other_dna.DNA[j-1]
                if not use_score_matrix:
                    if self_char == other_char:
                        diag_best = 1 + self.matching_matrix[i-1][j-1]
                else:
                    diag_best = self.score_matrix[self_char][other_char] + self.matching_matrix[i-1][j-1]
                left_best = self.matching_matrix[i-1][j] - penalty
                top_best = self.matching_matrix[i][j-1] - penalty
                best = max(diag_best, left_best, top_best)
                self.matching_matrix[i][j] = best
                if best == diag_best:
                    self.backtrack_matrix[i][j] = self.DIAGONAL
                elif best == left_best:
                    self.backtrack_matrix[i][j] = self.LEFT
                else:
                    self.backtrack_matrix[i][j] = self.TOP


    def graph_based_local_alignment(self, other_dna, penalty = 0, use_score_matrix = False):
        """
        Finds the best local alignment between the two sequences using dag best distance
        """
        adj_list = dict()
        self_len = len(self.DNA)
        other_len = len(other_dna.DNA)

        for i in range(self_len):
            for j in range(other_len):
                curr_node = (i, j)
                adj_list[curr_node] = [((i+1, j), -penalty), ((i, j+1), -penalty)]
                if not use_score_matrix:
                    if self.DNA[i] == other_dna.DNA[j]:
                        adj_list[curr_node] += [ ((i+1, j+1), 1) ]
                    else:
                        adj_list[curr_node] += [ ((i+1, j+1), 0) ]
                else:
                    adj_list[curr_node] += [ ((i+1, j+1), self.score_matrix[self.DNA[i]][other_dna.DNA[j]]) ]
            
        # Edges
        for i in range(self_len):
            adj_list[(i, other_len)] = [ ((i+1, other_len), -penalty) ]
        for i in range(other_len):
            adj_list[(self_len, i)] = [ ((self_len, i+1), -penalty) ]

        # adj_list is ready
        longest_path_neighbors = dict()
        longest_path_neighbors[(self_len+1, other_len+1)] = ((self_len+1, other_len+1), 0)
        for i in range(self_len+1, -1, -1):
            for j in range(other_len+1, -1, -1):
                curr_node = (i, j)
                if curr_node in longest_path_neighbors:
                    continue
                if curr_node not in adj_list:
                    longest_path_neighbors[curr_node] = (curr_node, 0)
                    continue
                best_neighbor = curr_node
                best_distance = 0
                for n in adj_list[curr_node]:
                    neighbor = n[0]
                    distance = n[1]
                    new_distance = longest_path_neighbors[neighbor][1]
                    if new_distance < 0:
                        continue
                    if new_distance + distance > best_distance:
                        best_distance = new_distance + distance
                        best_neighbor = neighbor
                longest_path_neighbors[(i, j)] = (best_neighbor, best_distance)


        best_starting_point = max(longest_path_neighbors, key=lambda x:longest_path_neighbors[x][1])
        best_score = longest_path_neighbors[best_starting_point][1]
        seq1 = ""
        seq2 = ""
        while best_starting_point != None and best_starting_point in longest_path_neighbors:
            next_best = longest_path_neighbors[best_starting_point][0]
            i = best_starting_point[0]
            j = best_starting_point[1]
            next_i = next_best[0]
            next_j = next_best[1]
            if i == next_i and j == next_j:
                break
            if i == next_i:
                seq1 += "-"
            else:
                seq1 += self.DNA[i]
            if j == next_j:
                seq2 += "-"
            else:
                seq2 += other_dna.DNA[j]
            best_starting_point = next_best
        
        return (best_score, seq1, seq2)
        

    def construct_subsequence(self, other_dna):
        """
        Constructs a subsequence from backtrack matrix
        """
        subseq = ""
        row = len(other_dna.DNA)
        col = len(self.DNA)
        while row > 0 and col > 0:
            if self.backtrack_matrix[col][row] == self.DIAGONAL:
                subseq = self.DNA[col-1] + subseq
                row -= 1
                col -= 1
            elif self.backtrack_matrix[col][row] == self.TOP:
                row -= 1
            else:
                col -= 1
        return subseq


    def construct_insdel_sequences(self, other_dna):
        """
        Replaces all ins-dels with '-' in matched sequences
        """
        seq1 = ""
        seq2 = ""
        row = len(other_dna.DNA)
        col = len(self.DNA)
        while row > 0 and col > 0:
            if self.backtrack_matrix[col][row] == self.DIAGONAL:
                seq1 = self.DNA[col-1] + seq1
                seq2 = other_dna.DNA[row-1] + seq2
                row -= 1
                col -= 1
            elif self.backtrack_matrix[col][row] == self.TOP:
                seq1 = "-" + seq1
                seq2 = other_dna.DNA[row-1] + seq2
                row -= 1
            else:
                seq1 = self.DNA[col-1] + seq1
                seq2 = "-" + seq2
                col -= 1

        if row == 0 and col != 0:
            while col > 0:
                seq1 = self.DNA[col-1] + seq1
                seq2 = "-" + seq2
		col -= 1
        elif row != 0 and col == 0:
            while row > 0:
                seq1 = "-" + seq1
                seq2 = other_dna.DNA[row-1] + seq2
		row -= 1
        return (seq1, seq2)
