from collections import defaultdict
from paired_kmer import PairedKmer

# a class to represent nodes in a debruijn graph
# slots:
#   _label: the paired k-1 mer label
#   _target_labels: list of target node labels for outgoing edges
class Node:
    # initialize node with label and empty targets list
    def __init__(self, label):
        self._label = label
        self._target_labels = []

    # return label for this node
    def label(self):
        return self._label

    # return list of target labels for this node
    def target_labels(self):
        return self._target_labels

    # get a string representation of node
    # "label -> <comma_separated_list_of_target_labels>"
    def __repr__(self):
        target_labels = [label.as_string() for label in self.target_labels()]
        targets_string = ",".join(target_labels)
        return self.label().as_string() + " -> " + targets_string

    # add target label to target list
    def add_target(self, target_label):
        self.target_labels().append(target_label)

    def remove_target(self, target_label):
        self.target_labels().remove(target_label)

    # return the number of targets for this node
    def num_targets(self):
        return len(self.target_labels())

# a class representing a Debruijn graph
# slots:
#   _nodes: a dictionary, keys are node labels,
#               values are objects of class Node
class Graph:
    # intialize with empty directory
    def __init__(self):
        self._nodes = {}

    # check if node with given label is in the node dictionary
    def __contains__(self, label):
        return label in self._nodes

    # return the node object corresponding to given label
    def __getitem__(self, label):
        return self._nodes[label]

    # iterate over nodes in graph
    def __iter__(self):
        return self._nodes.itervalues()

    # string representation of graph
    # calls __repr__ method for each node,
    # returns newline separated strings
    def __repr__(self):
        nodes = [node.__repr__() for node in self]
        return "\n".join(nodes)

    # add node to graph with given label
    # does not check if node with given label is already
    # in the graph, that has to be done elsewhere
    # creates a new object of class Node
    def add_node(self, label):
        self._nodes[label] = Node(label)

    # add an edge between nodes with given source and target labels
    # adds nodes to graph with corresponding label if needed
    def add_edge(self, source_label, target_label):
        if not source_label in self:
            self.add_node(source_label)
        if not target_label in self:
            self.add_node(target_label)
        self[source_label].add_target(target_label)

    def remove_edge(self, source_label, target_label):
        self[source_label].remove_target(target_label)

    # compute node degrees
    # return a dictionary with node labels as keys
    # values are tuples (in_degree, out_degree)
    def node_degrees(self):
        degrees = defaultdict(lambda: [0,0])
        for node in self:
            node_label = node.label()
            # set out-degree of node to the number of targets
            degrees[node_label][1] = node.num_targets()

            # increase the in-degree of targets by 1
            for target_label in node.target_labels():
                degrees[target_label][0] += 1
        return degrees

    # returns Nodes to precede given Node
    def get_ancestors(self, node):
        return [other_node for other_node in self if node.label() in other_node.target_labels()]

    # return Nodes preceeded by given Node
    def get_successors(self, node):
        return [self[target_label] for target_label in node.target_labels()]


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
        kmer = PairedKmer(kmer)
        source_label = kmer.prefix()
        target_label = kmer.suffix()

        # use the add edge method in graph class
        graph.add_edge(source_label, target_label)
    return graph
