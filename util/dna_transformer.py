"""
A class defining all transformations, like transcribing, base pairing etc. 
Remember that the DNA class member is read-only
"""
from dna_util import DNAUtil
from collections import defaultdict

class DNATransformer(DNAUtil):

    def __init__(self, input_dna):
        DNAUtil.__init__(self, input_dna)
        
    """
    Returns an RNA obtained after transcribing the DNA
    """
    def transcribe(self):
        ret = ""
        for c in self.DNA:
            if c == 'T':
                ret += "U"
            else:
                ret += str(c)
        return ret
            
    """
    Returns a complementary DNA strand for give DNA
    """
    def complement(self):
        ret = ""
        for c in self.DNA:
            ret += self.complement_dict[c]
        return ret
    
    """
    DNAs are usually stored as reverse complements. 
    This method finds the reverse complement of a DNA
    """
    def reverse_complement(self):
        ret = ""
        for c in self.DNA[::-1]:
            ret += self.complement_dict[c]
        return ret
