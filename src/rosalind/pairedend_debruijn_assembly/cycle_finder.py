from path import Path
from collections import defaultdict
import itertools

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
        self._target_iterators = dict([(node.label(), ProbedIter(node.target_labels())) for node in self._graph])

        # a map between nodes and items in paths being built
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
        it = self._target_iterators[source.label()]

        # check that there are more nodes
        if it.has_more():
            target_label = it.next()
            return self._graph[target_label]
        else:
            return None

    # check if node has unused out-going edges
    # input:
    #   node: graph node
    # returns:
    #   True or False
    def is_available(self, node):
        # out-going edge iterator for node
        it = self._target_iterators[node.label()]
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
