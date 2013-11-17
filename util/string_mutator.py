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

    
    def lengthKKmers(self, k):
        """
        Returns all length-k kmers
        """
        if (k == 1):
            return dna_util.bases
        ans = list()
        smallerKmers = self.lengthKKmers(k-1)
        for s in smallerKmers:
            for c in dna_util.bases:
                ans.append(c + s)
        return ans
            
