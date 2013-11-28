from util.dna_transformer import DNAMultiTransformer

class DNAGraphUtil(DNAMultiTransformer):

    def __init__(self, input_dnas, input_names=None):
        """
        Constructor keeps track of prefixes of the input dnas
        """
        DNAMultiTransformer.__init__(self, input_dnas, input_names)
        self.dna_prefixes = dict()
        for d in input_dnas:
            curr_prefix = d[:-1]
            try:
                self.dna_prefixes[curr_prefix][d] = 0
            except:
                self.dna_prefixes[curr_prefix] = dict()
		self.dna_prefixes[curr_prefix][d] = 0

 
    def adjacency_list(self):
        adj_list = dict()
        for d in self.dna_transformers:
            dna = d.DNA
            dna_suffix = dna[1:]
	    if dna_suffix in self.dna_prefixes:
                for other_dna in self.dna_prefixes[dna_suffix]:
                    try:
                        adj_list[dna] += [other_dna]
                    except:
                        adj_list[dna] = [other_dna]			    
        return adj_list
