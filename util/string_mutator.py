import copy
import dna_util

class StringMutator:
    
    def lengthKMutations(self, str, k):
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
