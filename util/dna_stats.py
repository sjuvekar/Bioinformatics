"""
A class containing all stats for input DNA, including the number of base-pairs, percents etc.
"""
import re
import dna_util
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
    Finds all occurrances of a motif in the DNA
    """
    def find_motif(self, text):
        return [m.start()+1 for m in re.finditer("(?={0})".format(text), self.DNA)]
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
    Returns the maximum GC contents among all DNAs
    """
    def max_gc_contents(self):
        max_contents = 0
        max_name = ""
        for dna in self.dna_stats:
            curr_contents = dna.gc_contents()
            if curr_contents > max_contents:
                max_contents = curr_contents
                max_name = dna.name
        return (max_contents, max_name)

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
    
    """
    Returns the consensus DNA, i.e. the DNA with bases occurring most frequently at every position
    """
    def consensus_dna(self):
        consensus_matrix = dict(zip(dna_util.bases, [[]] * len(dna_util.bases)))
        dna_length = len(self.dna_stats[0].DNA)
        for i in range(dna_length):
            curr_bases = defaultdict(int)
            for j in range(len(self.dna_stats)):
                curr_bases[self.dna_stats[j].DNA[i]] += 1
            for base in consensus_matrix.keys():
                consensus_matrix[base] = consensus_matrix[base] + [curr_bases[base]]
        consensus_string = ""
        for i in range(dna_length):
            consensus_string += max(consensus_matrix.keys(), key=lambda a: consensus_matrix[a][i])
        return (consensus_matrix, consensus_string)

    """
    Returns the longest_consecutive substring of all sequences
    """
    
