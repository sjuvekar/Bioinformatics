"""
A class defining all transformations, like transcribing, base pairing etc. 
Remember that the DNA class member is read-only
"""
import dna_util
from dna_util import DNAUtil
from collections import defaultdict

class DNATransformer(DNAUtil):

    def __init__(self, input_dna, fasta_name = "Rosalind_"):
        DNAUtil.__init__(self, input_dna, fasta_name)
        
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
            ret += dna_util.complement_dict[c]
        return ret
    
    """
    DNAs are usually stored as reverse complements. 
    This method finds the reverse complement of a DNA
    """
    def reverse_complement(self):
        ret = ""
        for c in self.DNA[::-1]:
            ret += dna_util.complement_dict[c]
        return ret


"""
A class to transform more than two DNAs, e.g. computing overlap graph
"""
class DNAMultiTransformer:

    """
    Constructor takes a list of DNAs and their names. 
    DNAs must all have same lengths, names can be empty
    """
    def __init__(self, input_dnas, input_names=None):
        if not input_names:
            input_names = [None] * len(input_dnas)
        self.dna_transformers = list()
        for i in range(len(input_dnas)):
            self.dna_transformers.append(DNATransformer(input_dnas[i], input_names[i]))

    
    """ 
    An O_k overlap graph, where there is an edge between (u, v) if k-suffix of DNA u matches k-prefix DNA v
    """
    def k_overlap(self, k):
        overlap_dict = dict()
        for i in range(len(self.dna_transformers)):
            for j in range(len(self.dna_transformers)):
                suff_dna = self.dna_transformers[i].DNA
                pref_dna = self.dna_transformers[j].DNA
                suff_name = self.dna_transformers[i].name
                pref_name = self.dna_transformers[j].name
                if suff_name == pref_name:
                    continue
                if suff_dna[-k:] == pref_dna[:k]:
                    try:
                        overlap_dict[suff_name] = overlap_dict[suff_name] + [pref_name]
                    except:
                        overlap_dict[suff_name] = [pref_name]
        return overlap_dict

