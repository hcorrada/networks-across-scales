import sys

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

        # init the edge hash for new node
        self._edges[nextnode] = {}
        
        return nextnode

    # splice a node between node and child with given label
    # split label at given length
    def splice(self, node, label, length):
        # find the node id for child with given label
        child = self.get_child(node, label)

        # get start end indices for edge label substring
        start, end = label

        # this is the position in text where the new node is spliced in
        splice = start + length - 1

        # add the spliced node as a child of node
        # with spliced label
        newlabel = (start, splice)
        newnode = self.add(node, newlabel)

        # now delete the existing edge
        del self._edges[node][label]

        # and add edge between spliced node and child
        newlabel = (splice+1, end)
        self._edges[newnode][newlabel] = child
        return newnode
        
    # process the ith suffix
    def process_suffix(self, suffix_index):
        # start at the trie's root
        curnode = self._root

        # we'll use i to walk through suffix starting
        # at suffix index
        i = suffix_index
        while i < len(self._text):
            # get matching edge label (and length of match)
            label, length = self.get_matching_edge(curnode, i)
            
            if label is None:
                # no matching edge, create a new child
                # label is [i,n] (n=len(text) - 1)
                newnode = self.add(curnode, (i, len(self._text)-1))

                # mark which suffix this leaf corresponds to
                self._leaf_labels[newnode] = suffix_index
                return

            # get (start,end) indices of labeling substring
            start, end = label

            # move pointer in suffix length positions
            i += length

            # if we used up all of the label in the current edge
            # then we set curnode to the corresponding child
            # if there is label left over we need to splice in a new node
            curnode = self.get_child(curnode, label) if (end - start + 1) == length else self.splice(curnode, label, length)

    def get_match_length(self, text):
        curnode = self._root

        i = 0
        total_length = 0
        
        while i < len(text):
            label, length = self.get_matching_edge(curnode, i, text)
            if label is None:
                return total_length

            total_length += length
            start, end = label
            i += length
            
            if (end - start + 1) == length:
                curnode = self.get_child(curnode, label)
            else:
                return total_length
        return total_length

    # find an outgoing edge of node with prefix of label matching a prefix of
    # text suffix at position 'start'
    #
    # if no matching prefix found, returns None,None
    # if matching prefix found, returns edgeLabel, matchLength
    # where matchLength is the length of matching prefix
    def get_matching_edge(self, node, start, text=None):
        if text is None:
            text = self._text
            
        # get the first character of suffix
        c = text[start]

        # get labels of all outgoing edges at node
        labels = self._edges[node].keys()

        # and the first character of each label
        chars = [self._text[x[0]] for x in labels]
        
        if c not in chars:
            # there is no matching first character,
            # i.e., no matching label
            return None, None

        # find the matching edge
        index = chars.index(c)
        label = labels[index]

        # now let's see how long is the matching prefix
        nlabel = label[1] - label[0] + 1
        length = 1
        while length < nlabel and start+length < len(text) and text[start+length] == self._text[label[0] + length]:
            length += 1

        # return the label of matching edge, and the length of the matching prefix
        return label, length
    
    # construct the tree
    def construct(self):
        # process each suffix in order
        for i in xrange(len(self._text)):
            self.process_suffix(i)
            
    # get the id of child of 'node' labeled by string label
    def get_child(self, node, label):
        # get the edge dictionary for node
        edges = self._edges[node]

        # return child id if there
        return edges[label] if label in edges else None

                  
# construct the suffix tree from a given text
def make_suffix_tree(text):
    # make an empty suffix tree
    tree = SuffixTree(text)

    # now process the suffixes
    tree.construct()
    return tree

def find_longest_common_substring(text, tree):
    best_index = -1
    best_length = 0
    for i in xrange(len(text)):
        cur_length = tree.get_match_length(text[i:])
        
        if cur_length >= best_length:
            best_length = cur_length
            best_index = i
    return text[best_index:(best_index+best_length)]

def readdat(filename):
    with open(filename, 'r') as f:
        text1 = f.readline().strip()
        text2 = f.readline().strip()
        return text1, text2

def main(filename):
    text1, text2 = readdat(filename)
    
    tree = make_suffix_tree(text1 + '$')
    lcs = find_longest_common_substring(text2, tree)
    print lcs
    
if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
    
