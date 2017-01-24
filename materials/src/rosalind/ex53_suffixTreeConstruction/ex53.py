import sys
from collections import deque

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

    # print the edge labels
    def print_labels(self):
        out = ''
        # iterate over each node
        for edges in self._edges.itervalues():
            # iterate over outgoing edges
            for label  in edges.iterkeys():
                out += (label + "\n")
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

    # compress non-branching edge labels
    def to_suffix_tree(self):
        # start the queue with root, empty label, and no parent
        q = deque([(self._root, '', None)])

        while len(q) > 0:
            # process the next node in queue
            # curlabel is the first character in conctenated string
            # parent is the first node in a potential chain
            curnode, curlabel, parent = q.popleft()

            # this variable keeps the concatenated string
            # we keep curlabel around to update the outgoing edge
            # on parent
            newlabel = curlabel

            # this is the current child in the potential chain
            child = curnode
            
            # get edges of node we just popped
            edges = self._edges[curnode]
            
            # while the current set of edges is non-branching
            while len(edges) == 1:
                if child != curnode:
                    # if this is an internal node in chain, remove it from trie
                    del self._edges[child]
                    
                # get the label and child for current edge set (of size 1)
                label, child = edges.items()[0]
                # append label to concatenated string
                newlabel += label

                # update the set of edges
                edges = self._edges[child]

            # check if we found a chain
            if child != curnode:
                # we did, remove the node we had popped from stack
                del self._edges[curnode]
                
                # update edges at the chain root (parent node)
                del self._edges[parent][curlabel]
                self._edges[parent][newlabel] = child

            # add outgoing edges at end of chain to queue if needed
            curnode = child
            if len(edges) > 1:
                    # update the q
                for label, child in edges.iteritems():
                    q.append((child, label, curnode))
         
         
# construct the trie from list of patterns
def make_suffix_trie(text):
    # make an empty trie
    trie = Trie()

    # insert each suffix
    for i in xrange(len(text)):
        trie.process(text[i:])
    return trie

def readdat(filename):
    with open(filename, 'r') as f:
        return f.readline().strip()

def main(filename):
    text = readdat(filename)
    trie = make_suffix_trie(text)
    trie.to_suffix_tree()
    print trie.print_labels()
    
if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
    
