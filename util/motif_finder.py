import copy
import numpy
import pandas
import dna_util
from dna_stats import DNAStats, DNAMultiStats
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
            
        
    def min_mismatch(self, kmer):
        """
        Finds the minimum mismatch berween kmger and any position in string
        """
        mismatch = len(self.DNA)
        kmer_len = len(kmer)
        for i in range(len(self.DNA)-kmer_len+1):
            mismatched_positions = self.mismatches(kmer, self.DNA[i:i+kmer_len])
            if mismatched_positions < mismatch:
                mismatch = mismatched_positions
        return mismatch


    def median_kmer(self, other_dnas, k):
        """
        Find the median kmer, i.e. the one reducing the distance from every DNA
        """
        best_distance = k * (1 + len(other_dnas))
        best_kmer = ""
        for kmer in self.stringMutator.lengthKKmers(k):
            distance = self.min_mismatch(kmer)
            for dna in other_dnas:
                distance += dna.min_mismatch(kmer)
            if distance < best_distance:
                best_distance = distance
                best_kmer = kmer
        return best_kmer


    def kmer_probability(self, kmer, profile):
        """
        Gives probability of kmer = product of probability of individual character. The profile is
        k X 4 pandas dataframe indexed by ['A', 'C', 'G', 'T'] columns
        """
        prob = 1
        for i in range(len(kmer)):
            char = kmer[i]
            prob *= profile[char][i]
        return prob

    
    def profile_most_probable_kmer(self, k, profile):
        """
        Find the kmer that has the maximum probability of appearing in profile. The profile is
        k X 4 pandas dataframe indexed by ['A', 'C', 'G', 'T'] columns
        """
        max_prob = 0
        best_kmer = self.DNA[0:k]
        for i in range(len(self.DNA)-k+1):
            kmer = self.DNA[i:i+k]
            prob = self.kmer_probability(kmer, profile)
            if prob > max_prob:
                max_prob = prob
                best_kmer = kmer
        return best_kmer


    def nucleotide_count(self, dnas, is_pseudo_count = False):
        """
        Computes a matrix of counts for every column in dnas
        """
        if is_pseudo_count:
            count_matrix = pandas.DataFrame(1., index = range(len(dnas[0])), columns = dna_util.bases)
        else:
            count_matrix = pandas.DataFrame(0., index = range(len(dnas[0])), columns = dna_util.bases)

        for i in range(len(dnas[0])):
            for j in range(len(dnas)):
                count_matrix[dnas[j][i]][i] += 1
        return count_matrix


    def profile(self, dnas, is_pseudo_count = False):
        """
        Creates a profile (prob. distribution) from dnas
        """
        count_matrix = self.nucleotide_count(dnas, is_pseudo_count)
        profile_matrix = count_matrix.copy()
        for i in profile_matrix.index:
            profile_matrix.ix[i] /= float(sum(profile_matrix.ix[i]))
        return profile_matrix

    
    def score(self, dnas):
        """
        Creates a hamming-distance score from consensus dna of given dnas
        """
        dna_multistat = DNAMultiStats(dnas)
        consensus_dna = DNAStats(dna_multistat.consensus_dna()[1])       
        score = 0.
        for d in dnas:
            score += consensus_dna.hamming_distance(DNAStats(d))
        return score

    
    def greedy_motif_search(self, dnas, k, is_pseudo_count = False):
        best_motif = [ self.DNA[0:k] ]
        for i in range(len(dnas)):
            best_motif.append(dnas[i].DNA[0:k])
        
        for i in range(len(self.DNA)-k+1):
            motifs = [self.DNA[i:i+k]]
            for j in range(len(dnas)):
                profile_so_far = self.profile(motifs, is_pseudo_count)
                new_motif = dnas[j].profile_most_probable_kmer(k, profile_so_far)
                motifs.append(new_motif)
            best_score = self.score(best_motif)
            new_score = self.score(motifs)
            if new_score < best_score:
                best_motif = copy.deepcopy(motifs)
        
        return best_motif
        

    def randomized_motif_search(self, dnas, k):
        i = numpy.random.randint(len(self.DNA)-k+1)
        best_motif = [ self.DNA[i:i+k] ]
        for i in range(len(dnas)):
            rand_i = numpy.random.randint(len(self.DNA)-k+1)
            best_motif.append(dnas[i].DNA[rand_i:rand_i+k])
        motif = copy.deepcopy(best_motif)
        while True:
            curr_profile = self.profile(motif, True)
            motif[0] = self.profile_most_probable_kmer(k, curr_profile)
            for i in range(len(dnas)):
                motif[i+1] = dnas[i].profile_most_probable_kmer(k, curr_profile)
            best_score = self.score(best_motif)
            new_score = self.score(motif)
            if new_score < best_score:
                best_motif = copy.deepcopy(motif)
            else:
                return best_motif
