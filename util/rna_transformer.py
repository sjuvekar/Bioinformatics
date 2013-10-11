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
        if len(self.RNA) <= 2:
            return ret
        count = 2
        while count < len(self.RNA):
            codon = self.RNA[count-2] + self.RNA[count-1] + self.RNA[count]
            if count == 2 and codon != rna_util.start_codon:
                return ""
            amino = rna_util.genetic_code[codon]
            if amino == "Stop":
                return ret
            ret += amino
            count += 3
        return ""
            
    """
    Find all possible translation of the current RNA. Any sequence that starts with a start codon, ends with stop codon
    and does not have any other intermediate stop codon is a candidate
    """
    def all_possible_translations(self):
        all_proteins = []
        for global_ptr in range(len(self.RNA)):
            protein = RNATransformer(self.RNA[global_ptr:]).translate()
            if protein != "":
                all_proteins.append(protein)
        return all_proteins
