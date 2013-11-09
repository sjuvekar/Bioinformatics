import protein_util
from protein_util import ProteinUtil
import copy
import heapq
from collections import defaultdict

class ProteinTransformer(ProteinUtil):
    """
    All transformations of a protein
    """

    def __init__(self, input_protein, fasta_name="Rosalind_"):
        ProteinUtil.__init__(self, input_protein, fasta_name)
        
    def possible_subpeptides(self):
        """
        Mass Spectrometer breaks the protein into its possible continuous subpeptides
        including empty and complete seq. Return a list of all such subpeptides. 
        ASSUMPTION: The protein is cyclic and the sub-peptides can "wrap around"
        """
        ret = [""]
        protein_len = len(self.protein)
        for l in range(1, protein_len):
            for i in range(protein_len):
                if i + l <= protein_len:
                    ret += [self.protein[i : i+l]]
                else:
                    ret += [self.protein[i:] + self.protein[:(i+l)%protein_len]]
        ret += [self.protein]
        return ret


    def cyclospectrum(self):
        """
        This is a integer mass spectrum of the protein, where masses are for every continuous 
        subpeptide of given protein. Assume that the protein is cyclic
        """
        subpeptides = self.possible_subpeptides()
        spectrum = map(lambda a: ProteinUtil(a).int_weight(), subpeptides)
        return sorted(spectrum)

    
    def cyclopeptide_sequence(self, spectrum):
        """
        Return all possible peptide sequences that are consistent with the given spectrum
        """
	# Make dict for faster access
	spectrum_dict = dict(zip(spectrum, spectrum))
        # Compute the max weight in spectrum. That has to be the weight of the acid
        max_spectrum_weight = max(spectrum)
        sequences = [""]
        final_sequences = []
        while len(sequences) > 0:
            # Copy sequences to new_seq
            new_seq = copy.deepcopy(sequences)
            sequences = list()
            # Add new amino acid and check for consistency
            for s in new_seq:
                sequence_weight = ProteinTransformer(s).int_weight()
                for a in protein_util.amino_acids_by_weight:
                    new_protein = s+a
                    new_weight = sequence_weight + protein_util.protein_mass_int[a]
                    # Case 1: Check if the new acid is consistent with spectrum
                    flag = True
                    for i in range(len(s)):
                        w = ProteinTransformer(s[i:]+a).int_weight()
                        if w not in spectrum_dict:
                            flag = False
                            break
                    if not flag:
                        continue
                    # Check if new weight attained maximum. If not, just append to pending sequences
                    elif new_weight < max_spectrum_weight:
                        sequences.append(new_protein)
                    # Finally, compare entire new spectrum with given spectrum
                    else:
                        new_spectrum = ProteinTransformer(new_protein).cyclospectrum()
                        if new_spectrum == spectrum:
                            final_sequences.append(new_protein)
                            
                                        
	# Done. Now return final sequences
        return final_sequences
            

    
    def leaderboard_cyclopeptide_sequence(self, leaderboard_len, spectrum, protein_masses):
        """
        This is a heuristic algorithm that keeps top-n approx matching spectrums to the given spectrum
        and returns the best scoring of all of them.
        """
        # Dict for fast access
        spectrum_dict = defaultdict(int)
        for s in spectrum:
            spectrum_dict[s] += 1

        # Keep the max spectrum weight
        max_spectrum_weight = max(spectrum)
        
        # Maintain the leaderboard
        leaderboard = [(0, list(), dict())]
        
        # Maintain the max score
        best_peptide = list()
        max_score = 0
        min_score = 0

        while len(leaderboard) > 0:
            old_max_score = max_score
            new_leaderboard = list()
            for old_seq in leaderboard:
                old_score = old_seq[0]
                old_acid = old_seq[1]
                old_acid_weight = sum(old_acid)
                for a in protein_masses:
                    old_dict = copy.deepcopy(old_seq[2])
                    new_acid = old_acid + [a]
                    new_acid_weight = old_acid_weight + a
                    if new_acid_weight > max_spectrum_weight:
                        continue
                    # Find intersection
                    suffixes = [old_acid[i:]+[a] for i in range(len(old_acid)+1)]
                    prefixes = [[a]+old_acid[:i] for i in range(len(old_acid)+1)]
                    prefixes.remove([a])
                    
                    for i in suffixes + prefixes:
                        a_weight = sum(i)
                        if a_weight in spectrum_dict:
                            if a_weight not in old_dict:
                                old_dict[a_weight] = 1
                            elif old_dict[a_weight] < spectrum_dict[a_weight]:
                                old_dict[a_weight] += 1
                                
                    score = sum(old_dict.values())
                    
                    if score >= max_score:
                        max_score = score
                        best_peptide = new_acid
                    if score >= min_score or len(new_leaderboard) < leaderboard_len:
                        heapq.heappush(new_leaderboard, (score, new_acid, old_dict))
                    if len(new_leaderboard) > leaderboard_len:
                        old_min = heapq.heappop(new_leaderboard)
                        min_score = min(new_leaderboard)[0]
                        if old_min[0] == max_score:
                            heapq.heappush(new_leaderboard, old_min)
                                                
            print max_score 
            leaderboard = copy.deepcopy(new_leaderboard)
            
        return best_peptide


    def convolution_spectrum(self, spectrum):
        """
        Computes the convolution spectrum of give spectrum
        """
        ret = []
        for i in range(len(spectrum)):
            for j in range(i+1, len(spectrum)):
                diff = abs(spectrum[j] - spectrum[i])
                if diff > 0:
                    ret.append(diff)
        return ret


    def convolution_cyclopeptide_sequence(self, spectrum, convolution_len, leaderboard_len):
        """
        Compute the sequence by 1) first computing convolution spectrum 2) taking convolution_len most
        frequent masses in spectrum and 3) Using leaderboard algorithm with only those masses
        """
        convolution_dict = defaultdict(int)
        for i in range(len(spectrum)):
            for j in range(i+1, len(spectrum)):
                diff = abs(spectrum[j] - spectrum[i])
                if diff >= 57 and diff <= 200:
                    convolution_dict[diff] += 1
   
        print convolution_dict
        # Sort the convolution spectrum find top-most entries
        sorted_convolution_dict = sorted(convolution_dict, key= lambda a: convolution_dict[a], reverse=True)
        convolution_list = list()
        for i in range(min(convolution_len, len(convolution_dict))):
            convolution_list.append(sorted_convolution_dict[i])
        # Add ties
        i = convolution_len
        while i < len(sorted_convolution_dict) and convolution_dict[sorted_convolution_dict[i]] == convolution_dict[sorted_convolution_dict[i-1]]:
            convolution_list.append(sorted_convolution_dict[i])
            i += 1
       
        # Done, now run leaderboard
        return self.leaderboard_cyclopeptide_sequence(leaderboard_len, spectrum, convolution_list)
