import sys
from collections import deque

# class representing the SuffixTree
class SuffixTree:
    def __init__(self, text):
        self._root = 1          # id of the root node
        self._lastid = 1        # last node id used
        self._edges = {}       # dictionary of edges, keys are node ids, values are dictionaries, keyed by edge label: (start,end) tuple
        self._edges[1] = {}   # edge dictionary for root
        self._text = text      # the text string
        self._leaf_labels = {}

    # print the tree
    def __str__(self):
        out = ""
        # iterate over each node
        for id, edges in self._edges.iteritems():
            # iterate over outgoing edges
            for (start,end), child in edges.iteritems():
                substr = self._text[start:(end+1)]
                out += "%d %d %s\n" % (id, child, substr)
        return out

    # print the edge labels
    def print_labels(self):
        out = ''
        # iterate over each node
        for edges in self._edges.itervalues():
            # iterate over outgoing edges
            for start, end  in edges.iterkeys():
                substr = self._text[start:(end+1)]
                out += (substr + "\n")
        return out

    # add a new node to the tree as a child of 'node'
    # with edge labeled by 'label'
    def add(self, node, label):
        # the id of the new node
        nextnode = self._lastid + 1

        # add the entry in the node's edge dictionary
        self._edges[node][label] = nextnode

        # increase the id generator
        self._lastid += 1

        self._edges[nextnode] = {}
        
        return nextnode

    # splice a node between node and child with given label
    # split label at given length
    def splice(self, node, label, length):
        child = self.get_child(node, label)

        start, end = label
        splice = start + length - 1
        newlabel = (start, splice)
        newnode = self.add(node, newlabel)
        del self._edges[node][label]

        newlabel = (splice+1, end)
        self._edges[newnode][newlabel] = child
        return newnode
        
    # process the ith suffix
    def process_suffix(self, suffix_index):
        # start at the trie's root
        curnode = self._root

        i = suffix_index
        while i < len(self._text):
            # get matching edge label (and length of match)
            label, length = self.get_matching_edge(curnode, i)
            
            if label is None:
                # no matching edge, create a new child
                newnode = self.add(curnode, (i, len(self._text)-1))
                leaf_labels[newnode] = suffix_index
                return

            start, end = label
            i += length
            curnode = self.get_child(curnode, label) if (end - start + 1) == length else self.splice(curnode, label, length)

    def get_matching_edge(self, node, start):
        c = self._text[start]
        labels = self._edges[node].keys()
        chars = [self._text[x[0]] for x in labels]
        
        if c not in chars:
            return None, None
        index = chars.index(c)
        label = labels[index]

        nlabel = label[1] - label[0] + 1
        length = 1
        while length < nlabel and start+length < len(self._text) and self._text[start+length] == self._text[label[0] + length]:
            length += 1
        return label, length
    
    # construct the tree
    def construct(self):
        for i in xrange(len(self._text)):
            self.process_suffix(i)
            
    # get the id of child of 'node' labeled by string label
    def get_child(self, node, label):
        # get the edge dictionary for node
        edges = self._edges[node]

        # return child id if there
        return edges[label] if label in edges else None

                  
# construct the trie from list of patterns
def make_suffix_tree(text):
    # make an empty trie
    tree = SuffixTree(text)
    tree.construct()
    return tree

def readdat(filename):
    with open(filename, 'r') as f:
        return f.readline().strip()

def main(filename):
    text = readdat(filename)
    tree = make_suffix_tree(text)
    print tree.print_labels()
    
if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
    
