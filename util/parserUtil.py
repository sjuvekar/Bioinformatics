"""
A parser for fasta Format. 
>Name
AA
AA
"""
class fastaParser:
    def __init__(self, inp):
        self.input = inp
        # Add a ficticious > at the end. This won't be used for parsing
        self.input.append(">")
        self.names = []
        self.sequences = []

    def parse(self):
        curr_name = ""
        curr_seq = ""
        for l in self.input:
            if l.startswith(">"):
                if curr_name != "" or curr_seq != "":
                    self.names.append(curr_name)
                    self.sequences.append(curr_seq)
                curr_name = l[1:]
                curr_seq = ""
            else:
                curr_seq += l
    
        
