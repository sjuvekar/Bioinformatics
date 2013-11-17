import dna_util
from oriC_finder import OriCFinder

class MotifFinder(OriCFinder):
    """
    A class designed to find the origin of replicator (oriC) in given DNA
    """
    def __init__(self, input_dna, fasta_name = "Rosalind_"):
        OriCFinder.__init__(self, input_dna, fasta_name)
        
    def k_d_motifs(self, k, d):
        """
        Returns a dict of all motifs of length k with at most d mutations from DNA
        """
        motifs = dict()
        for i in range(len(self.DNA)-k+1):
            kmer = self.DNA[i:i+k]
            for mutation in self.stringMutator.lengthKMutations(kmer, d):
                motifs[mutation] = 0
        return motifs


    def motif_enumeration(self, other_dnas, k, d):
        """
        Returns all motifs of length k with at most d mutations common in DNA and
        other_dnas
        """
        common_motifs = list()
        self_motif = self.k_d_motifs(k, d)
        other_motifs = list()
        for dna in other_dnas:
            other_motifs.append(dna.k_d_motifs(k, d))
        for m in self_motif:
            flag = True
            for other_m in other_motifs:
                if m not in other_m:
                    flag = False
                    break
            if flag:
                common_motifs.append(m)
        return common_motifs
            
        
