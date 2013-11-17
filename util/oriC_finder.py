import dna_util
from dna_transformer import DNATransformer
from string_mutator import StringMutator
from collections import defaultdict

class OriCFinder(DNATransformer):
    """
    A class designed to find the origin of replicator (oriC) in given DNA
    """
    def __init__(self, input_dna, fasta_name = "Rosalind_"):
        DNATransformer.__init__(self, input_dna, fasta_name)
        self.stringMutator = StringMutator()
        
    def frequent_kmers(self, k):
        """
        Finds the most frequent kemrs
        """
        d = defaultdict(int)
        for i in range(len(self.DNA)-k):
            key = self.DNA[i:i+k]
            d[key] += 1
        max_freq = max(d.values())
        kmers = list()
        for k in d.keys():
            if d[k] == max_freq:
                kmers.append(k)
        return kmers

    
    def pattern_positions(self, pattern):
        """
        Find all places in DNA where pattern occurs
        """
        positions = list()
        pat_len = len(pattern)
        for i in range(len(self.DNA)-pat_len):
            if self.DNA[i:i+pat_len] == pattern:
                positions.append(i)
        return positions


    def find_clumps(self, kmer_len, window_len, freq):
        """
        Find the clumps inside entire DNA having size window_len where there
        is atleast one seq of kmer_len occuring with freq
        """
        small_dna = OriCFinder(self.DNA[0:window_len])
        kmer_dict = small_dna.kmer_freq(kmer_len)
        freq_kmers_dict = dict()
        for k in kmer_dict.keys():
            if kmer_dict[k] >= freq:
                freq_kmers_dict[k] = 0
        # Now shift the window one character at a time and modify kmer_dict
        for i in range(window_len, len(self.DNA)):
            # expel old sequence
            expelled_ind = i - window_len
            expelled_seq = self.DNA[expelled_ind:expelled_ind+kmer_len]
            kmer_dict[expelled_seq] -= 1
            # Add new sequence
            new_seq = self.DNA[i-kmer_len+1:i+1]
            try:
                kmer_dict[new_seq] += 1
            except:
                kmer_dict[new_seq] = 1
            if kmer_dict[new_seq] >= freq:
                freq_kmers_dict[new_seq] = 0
        #Done. now return the dict
        return freq_kmers_dict
    

    def min_gc_skew(self):
        """
        Finds the positions with the minimum GC skew (#G - #C)
        """
        skew = 0
        minskew = len(self.DNA)
        skew_index_dict = dict()
        for i in range(len(self.DNA)):
            if self.DNA[i] == "C":
                skew -= 1
            elif self.DNA[i] == "G":
                skew +=1
            # Update min
            if skew < minskew:
                minskew = skew
            # update skew dict
            try:
                skew_index_dict[skew] += [i]
            except:
                skew_index_dict[skew] = [i]
        # Done. Now return the list of indices
        return skew_index_dict[minskew]


    def mismatches(self, str1, str2):
        ans = 0
        for i in range(len(str1)):
            if str1[i] != str2[i]:
                ans += 1
        return ans


    def approx_matches(self, pattern, mismatch_thresh):
        """
        Find approx matching positions of pattern in current DNA
        """

        matches = list()
        pattern_length = len(pattern)
        for i in range(len(self.DNA)-pattern_length+1):
            mismatched_positions = self.mismatches(pattern, self.DNA[i:i+pattern_length])
            if mismatched_positions <= mismatch_thresh:
                matches.append(i)
        return matches
                                              

    def frequent_words_with_mismatch(self, kmer_len, mismatch_thresh):
        """
        Find the most frequent kmers in this DNA with at most mismatch_thresh mutations
        """
        frequent_words = defaultdict(int)
        
        for i in range(len(self.DNA)-kmer_len+1):
            mutations = self.stringMutator.lengthKMutations(self.DNA[i:i+kmer_len], mismatch_thresh)
            for m in mutations:
                frequent_words[m] += 1
                        
        return frequent_words


    def frequent_words_with_mismatch_and_reverse_complement(self, kmer_len, mismatch_thresh):
        """
        Find the most frequent kmers in this DNA with at most mismatch_thresh mutations and counting reverse complements
        """
        frequent_words = defaultdict(int)
        
        for i in range(len(self.DNA)-kmer_len+1):
            mutations = self.stringMutator.lengthKMutations(self.DNA[i:i+kmer_len], mismatch_thresh)
            for m in mutations:
                frequent_words[m] += 1
                #Compute reverse complement
                reverse_complement = DNATransformer(m).reverse_complement()
                frequent_words[reverse_complement] += 1
                        
        return frequent_words
