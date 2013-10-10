"""
A class defining all transformations, like translating, mRNA etc. 
Remember that the RNA class member is read-only
"""
import rna_util
from rna_util import RNAUtil

class RNATransformer(RNAUtil):

    def __init__(self, input_dna, fasta_name = "Rosalind_"):
        RNAUtil.__init__(self, input_dna, fasta_name)
        
    """
    Returns a protein seqeunce obtained after translating the RNA
    """
    def translate(self):
        ret = ""
        count = 0
        while count < len(self.RNA):
            codon = self.RNA[count] + self.RNA[count+1] + self.RNA[count+2]
            amino = rna_util.genetic_code[codon]
            if amino != "Stop":
                ret += amino
            count += 3
        return ret
            
