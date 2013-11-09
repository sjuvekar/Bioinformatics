from collections import defaultdict
import rna_util

protein_mass = {
    "A":   71.03711,
    "C":   103.00919,
    "D":   115.02694,
    "E":   129.04259,
    "F":   147.06841,
    "G":   57.02146,
    "H":   137.05891,
    "I":   113.08406,
    "K":   128.09496,
    "L":   113.08406,
    "M":   131.04049,
    "N":   114.04293,
    "P":   97.05276,
    "Q":   128.05858,
    "R":   156.10111,
    "S":   87.03203,
    "T":   101.04768,
    "V":   99.06841,
    "W":   186.07931,
    "Y":   163.06333
}

protein_mass_int = {
    "A":   71,
    "C":   103,
    "D":   115,
    "E":   129,
    "F":   147,
    "G":   57,
    "H":   137,
    "I":   113,
    "K":   128,
    "L":   113,
    "M":   131,
    "N":   114,
    "P":   97,
    "Q":   128,
    "R":   156,
    "S":   87,
    "T":   101,
    "V":   99,
    "W":   186,
    "Y":   163
}

amino_acids_by_weight = ["G", "A", "S", "P", "V", "T", "C", "I", "L", "N", "D", "K", "Q", "E", "M", "H", "F", "R", "Y", "W"]

possible_rnas = defaultdict(int)
for v in rna_util.genetic_code.values():
    possible_rnas[v] += 1

class ProteinUtil(object):
    """
    This class defines all attributes related to an Proteins. Weights, folding etc
    """

    def __init__(self, input_protein, fasta_name="Rosalind_"):
        self.protein = input_protein
	self.name = fasta_name
        
    def weight(self):
        return sum(map(lambda a: protein_mass[a], self.protein))
    
    def int_weight(self):
        return sum(map(lambda a: protein_mass_int[a], self.protein))
    
    def possible_rna_sequences(self):
	ans = 1
	for aa in self.protein:
	    ans = (ans * possible_rnas[aa]) % 1000000
	ans = (ans * possible_rnas["Stop"]) % 1000000
	return ans
