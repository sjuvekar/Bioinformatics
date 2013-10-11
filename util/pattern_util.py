"""
A class for defining pattern matching data-structures: prefix tries, suffix trees etc
"""

"""
A prefix tries contains a dictionary {symbol -> prefix_trie} the maps a 
symbol on edge to a new trie node
"""
class PrefixTrie:
   
    def __init__(self, children_dict = None):
        self.id = 0
        self.children_dict = children_dict

    """
    Adds a sequence to the trie. 
    If there is an edge with symbol s starting the current node, then continue
    else, create an extra entry in the dict for new node corresponding to s and recurse
    """
    def add_sequence(self, seq):
        if seq == None or seq == "":
            return
        if not self.children_dict or seq[0] not in self.children_dict.keys():
            new_trie = PrefixTrie()
            new_trie.add_sequence(seq[1:])
            if not self.children_dict:
                self.children_dict = {}
            self.children_dict[seq[0]] = new_trie
        else:
            self.children_dict[seq[0]].add_sequence(seq[1:])

    
    """
    Give each node of the tree a unique id.
    Sets the id of current node to new_id and then start setting the ids of children from new_id + 1
    each child will return the id returned by children + 1
    """
    def set_id(self, new_id):
        self.id = new_id
        if not self.children_dict:
            return new_id + 1
        else:
            fresh_id = new_id + 1
            for child in self.children_dict.values():
                fresh_id = child.set_id(fresh_id)
            return fresh_id
        
    """
    Print: print (id -> child_id -> symbol)
    """
    def pprint(self):
        if not self.children_dict:
            return
        for (symbol, child) in self.children_dict.iteritems():
            print self.id, child.id, symbol
            child.pprint()
