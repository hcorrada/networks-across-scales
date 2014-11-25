import sys

# class representing the Trie
class Trie:
    def __init__(self):
        self._root = 1          # id of the root node
        self._lastid = 1        # last node id used
        self._edges = {}       # dictionary of edges, keys are node ids, values are dictionaries, keyed by edge label
        self._edges[1] = {}   # edge dictionary for root

    # print the trie
    def __str__(self):
        out = ""
        # iterate over each node
        for id, edges in self._edges.iteritems():
            # iterate over outgoing edges
            for label, child in edges.iteritems():
                out += "%d %d %s\n" % (id, child, label)
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
        return nextnode

    # insert a pattern to the trie
    def process(self, pattern):
        # turn the string into a list of characters
        charlist = list(pattern)

        # start at the trie's root
        curnode = self._root

        # for each character in pattern
        for c in charlist:
            # check if there is an outgoing edge
            # in current node labeled with character
            nextnode = self.get_child(curnode, c)

            if nextnode is None:
                # there isn't so add the outgoing edge
                nextnode = self.add(curnode, c)

                # add the node to the trie
                self._edges[nextnode] = {}
            # traverse to next node
            curnode = nextnode

    # get the id of child of 'node' labeled by character c
    def get_child(self, node, c):
        # get the edge dictionary for node
        edges = self._edges[node]

        # return child id if there
        return edges[c] if c in edges else None

# construct the trie from list of patterns
def make_trie(patterns):
    # make an empty trie
    trie = Trie()

    # insert each pattern
    for pattern in patterns:
        trie.process(pattern)
    return trie

def readdat(filename):
    with open(filename, 'r') as f:
        return [line.strip() for line in f]

def main(filename):
    patterns = readdat(filename)
    trie = make_trie(patterns)
    print trie
    
if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
    
