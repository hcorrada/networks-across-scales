from collections import defaultdict
import itertools
from collections import deque
from copy import deepcopy

# a class representing a node label
# using paired kmers
#
# slots:
#   _first: the first kmer in pair
#   _second: the second kmer in pair
class PairedKmer:
    # initialize from string with format
    # <first>|<second>
    def __init__(self, string, second=None):
        if second is None:
            self._first, self._second = string.split("|")
        else:
            self._first = string
            self._second = second

    def as_string(self):
        return self._first + "|" + self._second

    # return string representation of label with format
    # <first>|<second>
    def __repr__(self):
        return self.as_string()

    # how to hash this, we need this to
    # use this as dictionary key
    def __hash__(self):
        return hash(self.as_string())

    # equality based on string equality
    def __eq__(self, other):
        return self.as_string() == other.as_string()

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
#   _id: integer id assigned to node
#   _label: the node label, assumed to have a 'as_string' method
#   _targets: list of target nodes for outgoing edges
class Node:
    # initialize node with label, id and empty targets list
    def __init__(self, label, id):
        self._id = id
        self._label = label
        self._targets = []

    # return the id of this node
    def id(self):
        return self._id

    # two nodes are  equal if they have the same id
    def __eq__(self, other):
        return self.id() == other.id()

    # hash nodes based on id
    def __hash__(self):
        return hash(self.id())

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


    # get a string representation of node
    # "label -> <comma_separated_list_of_target_labels>"
    def __repr__(self):
        target_labels = [label.as_string() for label in self.target_labels()]
        targets_string = ",".join(target_labels)
        return self.label().as_string() + " -> " + targets_string

    # print using ids instead of labels
    def debug_print(self):
        me = str(self.id()) + ":" + self.label().as_string()
        return me + " -> " + self.target_ids().__repr__()

    # add target node to target list
    def add_target(self, target):
        self.targets().append(target)

    # remove target nodes from target list
    # based on target's id
    def remove_target(self, target):
        self.targets().remove(target)

    # return the number of targets for this node
    def num_targets(self):
        return len(self.targets())

# a class representing a graph
# slots:
#   _lastid: last id assigned to a node in the graph
#   _nodes: a dictionary of nodes, keys are id, values are objects of class Node
#   _label_map: a dictionary, keys are node labels,
#               values are node ids
class Graph:
    # intialize with empty directory
    def __init__(self):
        self._lastid = -1
        self._nodes = {}
        self._label_map = {}

    # check if node with given label is in the node dictionary
    def __contains__(self, label):
        return label in self._label_map

    # iterate over nodes in graph
    def __iter__(self):
        return self._nodes.itervalues()

    # string representation of graph
    # calls __repr__ method for each node,
    # returns newline separated strings
    def __repr__(self):
        nodes = [node.__repr__() for node in self]
        return "\n".join(nodes)

    # print using node ids instead of labels
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
        self._nodes[node.id()] = node

        if not label in self:
            self._label_map[label] = list()
        self._label_map[label].append(node.id())
        return node

    # get arbitrary node in graph with given label
    def get_node(self, label):
        node_id = self._label_map[label][0]
        return self._nodes[node_id]

    # get node in graph with given id
    def get_node_from_id(self, node_id):
        return self._nodes[node_id]

    # get arbitrary node in graph with given label,
    # if there is no node in graph with that label,
    # add a new node
    def get_or_make_node(self, label):
        return self.get_node(label) if label in self else self.add_node(label)

    # add an edge between given source and target nodes
    def add_edge(self, source, target):
        source.add_target(target)

    # remove edge between given source and target nodes
    def remove_edge(self, source, target):
        source.remove_target(target)

    # compute node degrees
    # return a dictionary with node ids as keys
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

    # return the total number of edges in graph
    def num_edges(self):
        # calculate degree for all nodes
        node_degrees = self.node_degrees()

        # add them up
        num_edges = 0
        for _, out_degree in node_degrees.itervalues():
            num_edges += out_degree
        return num_edges

    # returns Nodes that precede given Node
    def get_ancestors(self, node):
        return [other_node for other_node in self if node.id() in other_node.target_ids()]

    # return Nodes preceeded by given Node
    def get_successors(self, node):
        return node.targets()


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

# class representing item in a doubly-linked list
# we use this for path representation as a doubly-linked
# list of nodes
class DoubleList:
    def __init__(self, node, prev_item=None, next_item=None):
        self._node = node
        self._prev_item = prev_item
        self._next_item = next_item


# class representing a path
#
# slots:
#   _head: first item in doubly-linked list
#   _tail: last item in doubly_linked list
#   _item_map: mapping from node id to item in path
class Path:
    # set head and tail to None
    # and create an empty node->item map
    def __init__(self):
        self._head = None
        self._tail = None
        self._item_map = defaultdict(list)

    # return true if path is empty
    def is_empty(self):
        return self._head is None and self._tail is None

    # return the number of edges in the path
    def num_edges(self):
        # count the number of items
        num_nodes = 0
        for item_list in self._item_map.itervalues():
            num_nodes += len(item_list)
        # number of edges is the number of items minus 1
        return num_nodes - 1

    # append node to the tail of the list
    # add node to item map
    def append(self, node):
        # create a new item and add it
        # to the node->item map
        item = DoubleList(node)
        self._item_map[node.id()].append(item)

        if self.is_empty():
            # make new item head and tail of list
            self._head = self._tail = item
        else:
            # grab the item in the tail
            tail_item = self._tail

            # make item in tail point to new item
            tail_item._next_item = item

            # make new item's back pointer point to item in tail
            item._prev_item = tail_item

            # make the new item the tail of the list
            self._tail = item

    # find arbitrary item in path linked list containing given node
    def find_item(self, node):
        # grab the list of items for the node
        item_map = self._item_map[node.id()]
        if item_map is None or len(item_map) == 0:
            return None
        else:
            # return the first item in the list
            return item_map[0]

    # return list of items in path containing
    # given node
    def find_all_items(self, node):
        return self._item_map[node.id()]

    # extend this path with other path
    # at given item
    #
    # starting state:
    #   self: self.head <-> ... <-> prev <-> item <-> next ... <-> self.tail
    #   other: other.head <-> ... <-> other.tail
    #
    # after stitching we will have
    #   self.head <-> ... <-> prev <-> other.head <-> ... <-> other.tail <-> next <-> ... <-> self.tail
    def stitch(self, item, other):
        # if this path is empty, just
        # return the other path
        if self.is_empty():
            return other

        # if there is no place to stich the other path
        # just return this path
        if item is None:
            return self

        # get reference to item preceding stitch item on this path
        prev_item = item._prev_item

        # stick head of other path after preceding item on this path
        # state after this operation
        #   a) self.head <-> ... <-> prev <-> other.head <-> ... <-> other.tail
        #   b) item <-> next <-> ... <-> self.tail
        prev_item._next_item = other._head
        other._head._prev_item = prev_item

        # get reference to next item on this path
        next_item = item._next_item

        # put tail of other path before next item on this path
        # we get final state after this operation
        next_item._prev_item = other._tail
        other._tail._next_item = next_item

        # merge item maps
        for label, item_list in other._item_map.items():
            self._item_map[label] += item_list
        return self

    # a generator for nodes along the path
    def nodes(self):
        current = self._head
        while current is not None:
            yield current._node
            # make sure you stop after the tail
            # to avoid inifite loop
            if current == self._tail:
                current = None
            else:
                current = current._next_item

    # a string representation of the path
    def __repr__(self):
        node_labels = [node.label().as_string() for node in self.nodes()]
        return " -> ".join(node_labels)

# we use this helper method when reconstructing a string
# from a given path. Given a list of kmers, constructs a string
# by taking the first k-mer in the list and appending the last
# character of remaining k-mers in succession
#
# input:
#   kmers: a list of kmers
#
# output:
#   string spelled out by kmer list
def _string_helper(kmers):
    out = ""
    for kmer in kmers:
        # use the complete kmer if it's the first
        # i.e., when output string is empty
        if len(out) == 0:
            out = kmer
        else:
            # otherwise append last character of kmer
            out += kmer[-1]
    return out

# get string spelled out by path of paired k-mers
# output:
#   string spelled out by paired-kmers in path
def build_string(path, k, d):
    # get string spelled out by first k-mer in each pair
    first_kmers = [node.label().first() for node in path.nodes()]
    prefix_string = _string_helper(first_kmers)

    # get string spelled out by second k-mer in each pair
    second_kmers = [node.label().second() for node in path.nodes()]
    suffix_string = _string_helper(second_kmers)

    # check the overlap between prefix and suffix strings
    # is valid
    n = len(prefix_string)
    prefix_overlap = prefix_string[k + d:]
    suffix_overlap = suffix_string[:-(k + d)]

    if prefix_overlap != suffix_overlap:
        # invalid overlap, no string is spelled out by this path
        return None
    else:
        # valid overlap, return the string
        return prefix_string + suffix_string[-(k + d):]

# class implementing a probed iterator
# we use this to iterate over unused out-going edges in graph nodes
# and to check if there any unused edges left in the node
#
# it uses itertools.tee to create two iterators
# over out-going edges, one is one step ahead and is used
# to check for more available out-going edges
class ProbedIter:
    def __init__(self, it):
        self._it, self._next_it = itertools.tee(iter(it))
        self._advance()

    # move the look-ahead iterator forward
    def _advance(self):
        try:
            self._lookahead = self._next_it.next()
        except StopIteration:
            # no more items so set lookahead to None
            self._lookahead = None

    # return the next item in iterator
    def next(self):
        # move look ahead forward
        self._advance()

        # return next item in main iterator
        return self._it.next()

    # check if there are more items in iterator
    def has_more(self):
        # there is as long as lookahead is not None
        return self._lookahead is not None

# objects implementing Eulerian Cycle finding
# algorithm
#
# we use a class here to encapsulate all the auxiliary
# data structures used in the algorithm
class CycleFinder:
    def __init__(self, g):
        # store reference to graph
        self._graph = g

        # initialized iterators (with lookahead) over out-going edges
        # of each node
        self._target_iterators = dict([(node.id(), ProbedIter(node.targets())) for node in self._graph])

        # a map between node ids and items in paths being built
        self._cycle_items = defaultdict(list)

        # a set of nodes with unused out-going edges
        # that have been included in the cycle built so far
        self._available_nodes = set()

    # return the next target node from given node
    # input:
    #   source: node in graph
    # returns:
    #   target node in a unused out-going edge
    def next_target(self, source):
        # out-going edge iterator for source node
        it = self._target_iterators[source.id()]

        # check that there are more nodes
        if it.has_more():
            return it.next()
        else:
            return None

    # check if node has unused out-going edges
    # input:
    #   node: graph node
    # returns:
    #   True or False
    def is_available(self, node):
        # out-going edge iterator for node
        it = self._target_iterators[node.id()]
        # check if iterator has more out-going edges
        return it.has_more()

    # make a new cycle starting and ending at given node
    # input:
    #   start node: node with outgoing edges in graph
    def make_a_cycle(self, start_node):
        # start empty cycle
        cycle = Path()

        # set starting node as current node
        cur_node = start_node

        # keep going while current node has
        # an unused out-going edge
        while self.is_available(cur_node):
            cycle.append(cur_node)
            next_node = self.next_target(cur_node)

            # check current node has more unused out-going
            # edges
            if self.is_available(cur_node):
                self._available_nodes.add(cur_node)
            elif cur_node in self._available_nodes:
                # it has no more unused edges but is in
                # the list of available nodes, so remove it
                self._available_nodes.discard(cur_node)
            cur_node = next_node

        # we're in the same starting node so close the cycle
        cycle.append(cur_node)
        return cycle

    # find an Eulerian cycle in the graph
    def run(self):
        # initialize the cycle
        cycle = Path()

        # grab an arbitrary node to start
        # and add it to the set of nodes with
        # unused out-going edges
        node = iter(self._graph).next()
        self._available_nodes.add(node)

        # while there are any nodes with unused
        # out-going edges
        while len(self._available_nodes) > 0:
            # grab the next available node
            # and make a new cycle starting and ending
            # in that node
            node = self._available_nodes.pop()
            new_cycle = self.make_a_cycle(node)

            # find an item containing this node
            # in the current cycle linked list
            # and extend the with the newly found cycle
            # at that item
            item = cycle.find_item(node)
            cycle = cycle.stitch(item, new_cycle)
        return cycle

# a generator for simple graphs
# used as follows given an object graph of class Graph
#
# for simple_graph in SimpleGraphGenerator(graph):
#   do something with simple_graph
class SimpleGraphGenerator:
    def __init__(self, g):
        # the queue of current candidate graphs
        self._all_graphs = deque([g])

    # the generator function itself
    def __iter__(self):
        # run while there are any non-simple graphs left in queue
        while len(self._all_graphs) > 0:
            # get the next candidate graph in queue
            next_graph = deque.pop(self._all_graphs)

            # try to find a bypass node, will get None if there is
            # no node to bypass, i.e., the graph is simple
            node_to_bypass = self.find_bypass_node(next_graph)

            if node_to_bypass is None:
                # graph is simple, so yield it
                yield next_graph
            else:
                # graph is not simple, let's add all the bypass graphs
                # for the node to bypass to the queue
                for bypass_graph in self.get_bypass_graphs(next_graph, node_to_bypass):
                    # check connectedness later
                    self._all_graphs.append(bypass_graph)

    # find a node with in-degree greater than 1 to bypass
    def find_bypass_node(self, graph):
        # compute node degrees
        node_degrees = graph.node_degrees()

        # see if any node has in degree greater than one
        for node in graph:
            in_degree, _ = node_degrees[node.id()]
            if in_degree > 1:
                # here's one, return it
                return node
        # no bypass nodes found
        return None

    # generate all bypass graphs from given
    # bypass node
    def get_bypass_graphs(self, graph, node):
        # get all ancestors for the node
        # get all successors for the node
        ancestors = graph.get_ancestors(node)
        successors = graph.get_successors(node)

        # now, generate a bypass graph for each ancestor, successor pair
        for ancestor in ancestors:
            for successor in successors:
                # copy the graph
                new_graph = deepcopy(graph)

                # grab references to corresponding nodes in the new graph
                ancestor = new_graph.get_node_from_id(ancestor.id())
                successor = new_graph.get_node_from_id(successor.id())
                node = new_graph.get_node_from_id(node.id())

                # remove edges going through bypass node
                new_graph.remove_edge(ancestor, node)
                new_graph.remove_edge(node, successor)

                # add new node with same label as bypass node
                # and connect to ancestor and successor of bypass node
                new_node = new_graph.add_node(node.label())
                new_graph.add_edge(ancestor, new_node)
                new_graph.add_edge(new_node, successor)

                # done
                yield new_graph

# class implementing a generator of simple graphs
# used as follows given object g of class Graph
#
# for cycle in CycleGenerator(g):
#   do something with cycle
class CycleGenerator:
    def __init__(self, g):
        # use a simple graph generator
        self._it = SimpleGraphGenerator(g)
        # store the number of edges in graph, we'll use this to check
        # that cycles obtained from simple graphs use all edges
        # i.e., that simple graph is connected
        self._num_edges = g.num_edges()

    # the generator function itself
    def __iter__(self):
        # get simple graph from the simple graph generator
        for graph in self._it:
            # find Eulerian cycle in simple graph
            cycle = CycleFinder(graph).run()

            # check cycle uses all edges (i.e., simple graph is connected)
            if cycle.num_edges() == self._num_edges:
                # it does, yield the resulting cycle
                yield cycle

# find nodes corresponding to start and end of the path
# also checks for nearly balanced conditions to find
# an Eulerian path
#
# the start node must have in-degree = out-degree - 1
# the end node must have in-degree = out-degree + 1
# all other nodes must have in-degree = out-degree
#
# input:
#   g: object of class Graph
# returns:
#   tuple of nodes: start, end (if conditions are satisfied)
#   otherwise returns None
def _find_start_end_nodes(g):
    # compute node degrees
    node_degrees = g.node_degrees()

    # check balance conditions
    start = end = None
    for node in g:
        in_degree, out_degree = node_degrees[node.id()]

        # check if node is balanced
        if in_degree == out_degree:
            continue

        # check if this could be an end node and we
        # haven't found another end node yet
        if in_degree == out_degree + 1 and end is None:
            end = node
            continue

        # check if this could be a start node and we
        # haven't found another start node yet
        if in_degree == out_degree - 1 and start is None:
            start = node
            continue

        # node is unbalanced, doesn't qualify for start or end
        # or we already have start or end node
        # conditions not met
        return None
    # check we found both start and end nodes
    if start is None or end is None:
        return None
    return start, end

# remove added edge from end to start node
# to recover Eulerian path from Eulerian cycle
#
# input:
#   cycle: object of class Path
#   start: Node object for start node in path
#   end: Node object for end node in path
#
# returns:
#   object of class Path
def get_path_from_cycle(cycle, start, end):
    # the idea is to find the edge in the path
    # connecting end -> start

    # find all items in the path that correspond
    # to the end node
    end_items = cycle.find_all_items(end)

    # now let's find which of those connects
    # to the start node
    start_item = None
    end_item = None

    for item in end_items:
        end_item = item
        start_item = item._next_item

        if start_item._node == start:
            break
        else:
            start_item = None

    # make sure we did the right thing
    assert(start_item is not None)
    assert(end_item is not None)

    # first of all, this is a cycle
    # so the head and tail point to the same node
    # so let's splice out the copy pointed to by the
    # head pointer in the cycle
    cycle._tail._next_item = cycle._head._next_item

    # now start the path to return at the start item
    cycle._head = start_item

    # end the path at the end item
    cycle._tail = end_item

    # and make sure head and tail are valid
    cycle._head._prev_item = None
    cycle._tail._next_item = None

    return cycle

# preprocess graph before finding Eulerian cycle
#
# 1. checks balance conditions
# 2. finds start and end nodes
# 3. adds end->start edge
#
# input:
#   g: object of class Graph
# returns:
#   tuple if graph satisfies conditions: new_graph, start_node, end_node
#   None otherwise
def _preprocess_graph(g):
    # find start and end nodes
    # if balance conditions for graph are not satisfied
    # this will return None
    result = _find_start_end_nodes(g)

    # check if balancing conditions are met
    if result is None:
        return None

    # unpack start and end nodes
    start, end = result

    # add edge to connect start and end nodes
    g.add_edge(end, start)
    return g, start, end

# finds an Eulerian path, if it exists, in the given
# graph
#
# input:
#   g: object of class Graph
# returns:
#   object of class Path
def find_eulerian_path(g):
    # preprocess graph
    g, start, end = _preprocess_graph(g)

    # find Eulerian cycle
    cycle_finder = CycleFinder(g)
    cycle = cycle_finder.run()

    # remove extra edge from cycle to
    # recover Eulerican path
    path = get_path_from_cycle(cycle, start, end)

    # return path
    return path

# find valid eulerian path given graph
#
# input:
#   g: object of class Graph
#   validity_func: function Path->boolean that checks if generated path is valid
def find_valid_eulerian_path(g, validity_func):
    g, start, end = _preprocess_graph(g)

    # generate all Eulerian cycles
    for cycle in CycleGenerator(g):
        # get the next candidate cycle, turn into path
        path = get_path_from_cycle(cycle, start, end)
        if validity_func(path):
            # this is a valid path, so return it
            return path
    # couldn't find a valid path, so return None
    return None

import sys

# read input from file
# returns
#   paired_kmers: list of strings, each representing a pair of k-mers
#   k: k-mer length in pair
#   d: gap between k-mers in pair
def readdat(filename):
    with open(filename, 'r') as f:
        k, d = [int(s) for s in f.readline().strip().split()]
        paired_kmers = []
        for line in f:
            paired_kmers.append(line.strip())
        return paired_kmers, k, d


# main assembler routine for single paths
# prints assembled string, or returns None if
# single path does not generate a valid string
def single_path(filename):
    # read input from file
    paired_kmers, k, d = readdat(filename)

    # build DeBruijn graph from paired kmers
    g = build_graph(paired_kmers)

    # find Eulerian path in Debruijn graph
    # if there is no Eulerian path, this will return None
    path = find_eulerian_path(g)

    # check if an Eulerian path was found
    if path is None:
        return None

    # if so, construct a string from the given path
    # if path does not spell out a valid string
    # this will return None
    string = build_string(path, k, d)

    # print resulting string
    print string

# create a function that checks path validity
# it tries to construct a string from paired kmers
# and returns True if it can build a valid string
def make_validity_function(k, d):
    # try to build string from path
    # return True if successful
    def check_path(path):
        string = build_string(path, k, d)
        return string is not None
    return check_path

# main assembler routine where all possible paths
# may be generated, it generates paths until a valid is one
# is found, it returns None if no valid path is found
def all_paths(filename):
    # read input from file
    paired_kmers, k, d = readdat(filename)

    # build DeBruijn graph from paired kmers
    g = build_graph(paired_kmers)

    # find Eulerian path that generates valid string
    # if there is no valid Eulerian path, this will return None
    path = find_valid_eulerian_path(g, make_validity_function(k,d))

    # check if a valid path was found
    if path is None:
        return None

    # get string from valid path
    string = build_string(path, k, d)
    print string

# this is here to play nicely with ipython
if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    all_paths(filename)
