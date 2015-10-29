from collections import deque
from copy import deepcopy
from cycle_finder import CycleFinder

class SimpleGraphGenerator:
    def __init__(self, g):
        # the queue of current candidate graphs
        self._all_graphs = deque([g])

    def __iter__(self):
        while len(self._all_graphs) > 0:
            next_graph = deque.pop(self._all_graphs)
            node_to_bypass = self.find_bypass_node(next_graph)
#            print node_to_bypass

            if node_to_bypass is None:
                yield next_graph
            else:
                for bypass_graph in self.get_bypass_graphs(next_graph, node_to_bypass):
                    # check connectedness later
                    self._all_graphs.append(bypass_graph)

    def is_connected(self, graph):
        return True

    def find_bypass_node(self, graph):
        node_degrees = graph.node_degrees()

        for node in graph:
            in_degree, _ = node_degrees[node.id()]
            if in_degree > 1:
                return node
        return None

    def get_bypass_graphs(self, graph, node):
        ancestors = graph.get_ancestors(node)
        successors = graph.get_successors(node)

        for ancestor in ancestors:
            for successor in successors:
#                print "making bypass"
#                print node.debug_print()
#                print ancestor.debug_print()
#                print successor.debug_print()

                new_graph = deepcopy(graph)
#                print "graph copy"
#                print new_graph.debug_print()

                ancestor = new_graph.get_node_from_id(ancestor.id())
                successor = new_graph.get_node_from_id(successor.id())
                node = new_graph.get_node_from_id(node.id())

                new_graph.remove_edge(ancestor, node)
                new_graph.remove_edge(node, successor)

#                print "new graph"
#                print new_graph.debug_print()

                new_node = new_graph.add_node(node.label())
#                print "new node"
#                print new_node.debug_print()

                new_graph.add_edge(ancestor, new_node)
                new_graph.add_edge(new_node, successor)
#                print new_graph.debug_print()

                yield new_graph

class CycleGenerator:
    def __init__(self, g):
        self._it = SimpleGraphGenerator(g)
        self._num_edges = g.num_edges()

    def __iter__(self):
        for graph in self._it:
            cycle = CycleFinder(graph).run()
            # check cycle uses all edges
#            print cycle
            if cycle.num_edges() == self._num_edges:
                yield cycle
