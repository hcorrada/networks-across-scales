from collections import deque
from copy import deepcopy
from cycle_finder import CycleFinder

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
