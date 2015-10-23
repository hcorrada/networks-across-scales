# a class representing a node label
# using paired kmers
#
# slots:
#   _first: the first kmer in pair
#   _second: the second kmer in pair
class PairedKmer:
    # initialize from string with format
    # <first>|<second>
    def __init__(self, string):
        self._first, self._second = string.split("|")

    # return string representation of label with format
    # <first>|<second>
    def __repr__(self):
        return self._first + "|" + self._second

    # return the first kmer in pair
    def first(self):
        return self._first

    # return the second kmer in pair
    def second(self):
        return self._second

    # return the paired k-1 mer prefix
    def prefix(self):
        return PairedKmer(self._first[:-1], self._second[:-1])

    # return the paired k-1 mer suffix
    def suffix(self):
        return PairedKmer(self._first[1:], self._second[1:])

# a class to represent nodes in a debruijn graph
# slots:
#   _label: the k-1 mer label
#   _targets: list of target nodes for outgoing edges
class Node:
    # initialize node with label and empty targets list
    def __init__(self, label):
        self._label = label
        self._targets = []

    def label(self):
        return self._label

    def label_string(self):
        return self._label.__repr__

    def targets(self):
        return self._targets

    # get a string representation of node
    # "label -> <comma_separated_list_of_target_labels>"
    def __repr__(self):
        target_labels = [target.label_string() for target in self.targets()]
        targets_string = ",".join(targets)
        return self._label + " -> " + targets_string

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

# build paired debruijn graph, i.e., k,d-mer overlap graph
# from given set of kmers
#
# input:
#   kmers: list of kmers
# output:
#   an object of class Graph
def build_graph(paired_kmers):
    # initalize graph object
    graph = Graph()

    # add an edge for each kmer in list
    for kmer in paired_kmers:
        source_label = kmer.prefix()
        target_label = kmer.suffix()

        # use the add edge method in graph class
        graph.add_edge(source_label, target_label)
    return graph
