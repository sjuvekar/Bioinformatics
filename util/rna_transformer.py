import rna_util
from rna_util import RNAUtil

class RNATransformer(RNAUtil):
    """
    A class defining all transformations, like translating, mRNA etc. 
    Remember that the RNA class member is read-only
    """

    def __init__(self, input_dna, fasta_name = "Rosalind_"):
        RNAUtil.__init__(self, input_dna, fasta_name)


    def reverse_transcribe(self):
        """
        Reverse transcribe the RNA molecule to DNA
        """
        ans = ""
        for r in self.RNA:
            if r == "U":
                ans += "T"
            else:
                ans += r
        return ans


    def reverse_complement(self):
        """
        Returns the reverse complement of the RNA
        """
        ret = ""
        for r in self.RNA[::-1]:
            ret += rna_util.complement_dict[r]
        return ret


    def translate(self):
        """
        Returns a protein seqeunce obtained after translating the RNA
        """
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
            

    def all_possible_translations(self):
        """
        Find all possible translation of the current RNA. Any sequence that starts with a start codon, ends with stop codon
        and does not have any other intermediate stop codon is a candidate
        """
        all_proteins = []
        for global_ptr in range(len(self.RNA)):
            protein = RNATransformer(self.RNA[global_ptr:]).translate()
            if protein != "":
                all_proteins.append(protein)
        return all_proteins


    def does_encode(self, amino_acid):
        """
        Returns true if the RNA encodes for the given amino acid
        """
        count = 2
        acid_count = 0
        while count < len(self.RNA):
            codon = self.RNA[count-2] + self.RNA[count-1] + self.RNA[count]
            if rna_util.genetic_code[codon] != amino_acid[acid_count]:
                return False
            acid_count += 1
            count += 3
        return True

    def encoding_strings(self, amino_acid):
        """
        Returns all amino_acid encoding substrings of RNA
        """
        encoders = list()
        amino_length = len(amino_acid)
        for i in range(len(self.RNA) - 3 * amino_length + 1):
            # Check for forward strand
            curr_rna = RNATransformer(self.RNA[i : i + 3 * amino_length])
            if curr_rna.does_encode(amino_acid):
                encoders.append(curr_rna.RNA)
            # Check for reverse strand
            rev_rna = RNATransformer(curr_rna.reverse_complement())
            if rev_rna.does_encode(amino_acid):
                encoders.append(curr_rna.RNA)
        return encoders
