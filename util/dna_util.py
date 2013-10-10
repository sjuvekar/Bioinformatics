bases = ['A', 'C', 'G', 'T']
complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        
class DNAUtil(object):
    
    def __init__(self, input_dna, fasta_name="Rosalind_"):
        self.DNA = input_dna
	self.name = fasta_name
        
