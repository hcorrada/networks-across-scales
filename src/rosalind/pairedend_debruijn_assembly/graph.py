from collections import defaultdict
from paired_kmer import PairedKmer

# a class to represent nodes in a debruijn graph
# slots:
#   _label: the paired k-1 mer label
#   _target_labels: list of target node labels for outgoing edges
class Node:
    # initialize node with label and empty targets list
    def __init__(self, label, id):
        self._id = id
        self._label = label
        self._targets = []

    def __eq__(self, other):
        return self._id == other._id

    def __hash__(self):
        return hash(self._id)

    # return label for this node
    def label(self):
        return self._label

    # return list of target Node objects
    def targets(self):
        return self._targets

    # return list of target ids for this node
    def target_ids(self):
        return [target.id() for target in self.targets()]

    # return list of target labels
    def target_labels(self):
        return [target.label() for target in self.targets()]

    def id(self):
        return self._id

    # get a string representation of node
    # "label -> <comma_separated_list_of_target_labels>"
    def __repr__(self):
        target_labels = [label.as_string() for label in self.target_labels()]
        targets_string = ",".join(target_labels)
        return self.label().as_string() + " -> " + targets_string

    def debug_print(self):
        me = str(self.id()) + ":" + self.label().as_string()
        target_ids = [str(id) for id in self.target_ids()]
        return me + " -> " + target_ids.__repr__()

    # add target label to target list
    def add_target(self, target):
        self.targets().append(target)

    def remove_target(self, target):
        self.target_ids().remove(target)

    # return the number of targets for this node
    def num_targets(self):
        return len(self.targets())

# a class representing a Debruijn graph
# slots:
#   _nodes: a dictionary, keys are node labels,
#               values are objects of class Node
class Graph:
    # intialize with empty directory
    def __init__(self):
        self._lastid = -1
        self._nodes = []
        self._label_map = {}

    # check if node with given label is in the node dictionary
    def __contains__(self, label):
        return label in self._label_map

    # iterate over nodes in graph
    def __iter__(self):
        return iter(self._nodes)

    # string representation of graph
    # calls __repr__ method for each node,
    # returns newline separated strings
    def __repr__(self):
        nodes = [node.__repr__() for node in self]
        return "\n".join(nodes)

    def debug_print(self):
        nodes = [node.debug_print() for node in self]
        return "\n".join(nodes)

    # add node to graph with given label
    # does not check if node with given label is already
    # in the graph, that has to be done elsewhere
    # creates a new object of class Node
    def add_node(self, label):
        self._lastid += 1
        node = Node(label, self._lastid)
        self._nodes.append(node)

        if not label in self:
            self._label_map[label] = list()
        self._label_map[label].append(node.id())
        return node

    def get_node(self, label):
        node_id = self._label_map[label][0]
        return self._nodes[node_id]

    def get_or_make_node(self, label):
        return self.get_node(label) if label in self else self.add_node(label)

    # add an edge between given source and target nodes
    def add_edge(self, source, target):
        source.add_target(target)

    def remove_edge(self, source, target):
        source.remove_target(target)

    # compute node degrees
    # return a dictionary with node labels as keys
    # values are tuples (in_degree, out_degree)
    def node_degrees(self):
        degrees = defaultdict(lambda: [0,0])
        for node in self:
            node_id = node.id()
            # set out-degree of node to the number of targets
            degrees[node_id][1] = node.num_targets()

            # increase the in-degree of targets by 1
            for target_id in node.target_ids():
                degrees[target_id][0] += 1
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

        # grab source and target nodes if they exist,
        # if not create new nodes
        source = graph.get_or_make_node(source_label)
        target = graph.get_or_make_node(target_label)


        # use the add edge method in graph class
        graph.add_edge(source, target)
    return graph
