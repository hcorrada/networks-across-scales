from collections import deque
from copy import deepcopy

class SimpleGraphGenerator:
    def __init__(self, g):
        # the queue of current candidate graphs
        self._all_graphs = deque([g])

    def __iter__(self):
        while len(self._all_graphs) > 0:
            next_graph = deque.pop(self._all_graphs)
            node_to_bypass = self.find_bypass_node(next_graph)
            print node_to_bypass

            if node_to_bypass is None:
                yield next_graph
            else:
                for bypass_graph in self.get_bypass_graphs(next_graph, node_to_bypass):
                    if self.is_connected(bypass_graph):
                        self._all_graphs.append(bypass_graph)

    def find_bypass_node(self, graph):
        node_degrees = graph.node_degrees()

        for node in graph:
            in_degree, _ = node_degrees[node.label()]
            if in_degree > 1:
                return node
        return None

    def get_bypass_graphs(self, graph, node):
        ancestors = graph.get_ancestors(node)
        successors = graph.get_successors(node)

        for ancestor in ancestors:
            for successor in successors:
                print node
                print ancestor
                print successor

                new_graph = deepcopy(graph)
                new_graph.remove_edge(ancestor, node)
                new_graph.remove_edge(node, successor)

                print new_graph

                new_node = Node(node.label())
                new_graph.add_edge(ancestor, new_node)
                new_graph.add_edge(new_node, successor)
                yield new_graph

class CycleGenerator:
    def __init__(self, g):
        self._it = SimpleGraphGenerator(g)

    def __iter__(self):
        for graph in self._it:
            yield CycleFinder(graph).run()
