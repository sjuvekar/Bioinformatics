import bisect
import copy
import dna_util

class StringMutator:
    
    def lengthKMutations(self, str, k):
        """
        Mutations of a give strings with all strings differing from str by at most k letters
        """
        if k == 0 or str =="":
            return [str]
        mutations = []
        
        # First keep first character and recurse over smaller string
        smallerStringMutations = self.lengthKMutations(str[1:], k)
        for s in smallerStringMutations:
            mutations.append(str[0] + s)

        # Next mutate first letter and recurse
        new_bases = copy.deepcopy(dna_util.bases)
        new_bases.remove(str[0])
        smallerStringMutations = self.lengthKMutations(str[1:], k-1)
        for s in smallerStringMutations:
            for n in new_bases:
                mutations.append(n + s)
        
        # Return mutations
        return mutations

    
    def lengthKKmers(self, k, bases=dna_util.bases):
        """
        Returns all length-k kmers
        """
        if (k == 1):
            return bases
        ans = list()
        smallerKmers = self.lengthKKmers(k-1, bases)
        for s in smallerKmers:
            for c in bases:
                ans.append(c + s)
        return ans
            

    def lexicographic_kmers(self, string, k):
        """
        Returns all kmers of length k lexicographically sorted
        """
        ret = list()
        for i in range(len(string)-k+1):
            kmer = string[i:i+k]
            index = bisect.bisect_left(ret, kmer)
            ret.insert(index, kmer)
        return ret
