"""
A class containing all stats for input DNA, including the number of base-pairs, percents etc.
"""
from dna_util import DNAUtil
from collections import defaultdict

class DNAStats(DNAUtil):

    def __init__(self, input_dna, fasta_name="Rosalind_"):
        DNAUtil.__init__(self, input_dna, fasta_name)

    """
    Returns a dict containing counts of all nucleotides
    """
    def count_bases(self):
        ret = defaultdict(int)
        for c in self.DNA:
            ret[c] += 1
        return ret
            
    """ 
    Returns of percent of GC contents in the DNA
    """
    def gc_contents(self):
        d = self.count_bases()
        num = d['C'] + d['G']
        den = sum(d.values())
        return float(num) / float(den)

"""
A class containing SNA stats over multiple DNAs, like hamming distance, alignment etc
"""
class DNAMultiStats:
    
    """
    Constructor takes a list of DNAs and their names. 
    DNAs must all have same lengths, names can be empty
    """
    def __init__(self, input_dnas, input_names=None):
        if not input_names:
            input_names = [None] * len(input_dnas)
        self.dna_stats = list()
        for i in range(len(input_dnas)):
            self.dna_stats.append(DNAStats(input_dnas[i], input_names[i]))

    
    """
    Returns number of places where all the sequences differ
    """
    def hamming_distance(self):
       n_mismatch = 0
       dna_length = len(self.dna_stats[0].DNA)
       for i in range(dna_length):
           curr_base = self.dna_stats[0].DNA[i]
           for j in range(1, len(self.dna_stats)):
               if self.dna_stats[j].DNA[i] != curr_base:
                   n_mismatch += 1
                   break
       return n_mismatch
    

