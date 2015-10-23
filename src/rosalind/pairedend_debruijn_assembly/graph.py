import sys

# compute the kmer composition of string
# input:
#   k: k-mer length
#   text: dna string
# output:
#   list of k-mers in text
def kmer_composition(k, text):
    kmers = []
    for i in xrange(len(text) - k +1):
        kmers.append(text[slice(i,i+k)])
    kmers.sort()
    return kmers

# read input from file
# expected format:
#   k
#   text
# output:
#   tuple k, text
def readdat(filename):
    with open(filename, 'r') as f:
        k = int(f.readline().strip())
        test = f.readline().strip()
    return k, test

# a class to represent nodes in a debruijn graph
# slots:
#   _label: the k-1 mer label
#   _targets: list of target nodes for outgoing edges
class Node:
    # initialize node with label and empty targets list
    def __init__(self, label):
        self._label = label
        self._targets = []

    # get a string representation of node
    # "label -> <comma_separated_list_of_target_labels>"
    def __repr__(self):
        targets = [target._label for target in self._targets]
        targets = ",".join(targets)
        return self._label + " -> " + targets

    # add target node to target list
    def add_target(self, target):
        self._targets.append(target)

    # return the number of targets for this node
    def num_targets(self):
        return len(self._targets)

# a class representing a Debruijn graph
# slots:
#   _nodes: a dictionary, keys are node labels, values are Node objects
class Graph:
    # intialize with empty directory
    def __init__(self):
        self._nodes = {}

    # check if node with given label is in the node dictionary
    def __contains__(self, label):
        return label in self._nodes

    # return the node corresponding to given label
    def __getitem__(self, label):
        return self._nodes[label]

    # iterate over nodes in graph
    def __iter__(self):
        return self._nodes.itervalues()

    # string representation of graph
    # calls __repr__ method for each node,
    # sorts resulting strings
    # returns newline separated strings
    def __repr__(self):
        nodes = [node.__repr__() for node in self if node.num_targets() > 0]
        return "\n".join(sorted(nodes))

    # add node to graph with given label
    # does not check if node with given label is already
    # in the graph
    def add_node(self, label):
        self._nodes[label] = Node(label)

    # add an edge between nodes with given source and target labels
    # adds nodes to graph with corresponding label if needed
    def add_edge(self, source_label, target_label):
        if not source_label in self:
            self.add_node(source_label)
        if not target_label in self:
            self.add_node(target_label)
        self[source_label].add_target(self[target_label])

# build debruijn graph, i.e., k-mer overlap graph
# from given set of kmers
#
# input:
#   kmers: list of kmers
# output:
#   an object of class Graph
def build_graph(kmers):
    # initalize graph object
    graph = Graph()

    # add an edge for each kmer in list
    for kmer in kmers:
        source_label = kmer[:-1]
        target_label = kmer[1:]

        # use the add edge method in graph class
        graph.add_edge(source_label, target_label)
    return graph

def main(filename):
    k, text = readdat(filename)
    kmers = kmer_composition(k, text)
    graph = build_graph(kmers)
    print graph

# this is here so this plays nicely with ipython %loadpy magic
if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
