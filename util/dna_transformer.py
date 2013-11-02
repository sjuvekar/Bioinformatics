import dna_util
from dna_util import DNAUtil
from collections import defaultdict

class DNATransformer(DNAUtil):
    """
    A class defining all transformations, like transcribing, base pairing etc. 
    Remember that the DNA class member is read-only
    """

    def __init__(self, input_dna, fasta_name = "Rosalind_"):
        DNAUtil.__init__(self, input_dna, fasta_name)
        
    def transcribe(self):
        """
        Returns an RNA obtained after transcribing the DNA
        """
        ret = ""
        for c in self.DNA:
            if c == 'T':
                ret += "U"
            else:
                ret += str(c)
        return ret
            
    def complement(self):
        """
        Returns a complementary DNA strand for give DNA
        """
        ret = ""
        for c in self.DNA:
            ret += dna_util.complement_dict[c]
        return ret
    
    def reverse_complement(self):
        """
        DNAs are usually stored as reverse complements. 
        This method finds the reverse complement of a DNA
        """
        ret = ""
        for c in self.DNA[::-1]:
            ret += dna_util.complement_dict[c]
        return ret

    def splice(self, introns):
        """
        returns a spliced DNA strand obtained after splicing a list of interons from current DNA
        Assume only one such string exists
        """
        ret = ""
        ind = 0
        while ind < len(self.DNA):
            flag = False
            for i in introns:
                if self.DNA[ind:].startswith(i):
                    ind += len(i)
                    flag = True
                    break
            if not flag:
                ret += str(self.DNA[ind])
                ind += 1
        return ret
                

class DNAMultiTransformer:
    """
    A class to transform more than two DNAs, e.g. computing overlap graph
    """

    def __init__(self, input_dnas, input_names=None):
        """
        Constructor takes a list of DNAs and their names. 
        DNAs must all have same lengths, names can be empty
        """
        if not input_names:
            input_names = [None] * len(input_dnas)
        self.dna_transformers = list()
        for i in range(len(input_dnas)):
            self.dna_transformers.append(DNATransformer(input_dnas[i], input_names[i]))

    
    def k_overlap(self, k):
        """ 
        An O_k overlap graph, where there is an edge between (u, v) if k-suffix of DNA u matches k-prefix DNA v
        """
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

