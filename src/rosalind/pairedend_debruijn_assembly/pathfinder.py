from path import Path
from collections import defaultdict

class CycleFinder:
    def __init__(self, g, k, d):
        self._graph = g
        self._k = k
        self._d = d

        self._available = dict([(node.label(), True) for node in self._graph])
        self._target_iterators = dict([(node.label(), iter(node.target_labels())) for node in self._graph])
        self._cycle_items = defaultdict(list)

    def make_a_cycle(self, start_node):
        path = Path(self._k, self._d)
        path.append(start_node)
        return path

    def find_cycle(self):
        cycle = Path(self._k, self._d)

        for node in self._graph:
            if not self._available[node.label()]:
                continue

            self._available[node.label()] = False

            print "next node:"
            print node
            
            new_cycle = self.make_a_cycle(node)
            print new_cycle
            cycle = cycle.stitch(None, new_cycle)
            print cycle
            print

        return cycle

def _find_start_end_nodes(g):
    # compute node degrees
    node_degrees = g.node_degrees()

    # check balance conditions
    start = end = None
    for node in g:
        in_degree, out_degree = node_degrees[node.label()]

        if in_degree == out_degree:
            continue

        if in_degree == out_degree + 1 and end is None:
            end = node
            continue

        if in_degree == out_degree - 1 and start is None:
            start = node
            continue

        return None
    if start is None or end is None:
        return None
    return start, end

def find_eulerian_path(g, k, d):
    # find start and end nodes (also checks balance conditions)
    result = _find_start_end_nodes(g)
    if result is None:
        return None

    start, end = result

    # add edge to connect start and end node
    g.add_edge(end.label(), start.label())

    # get eulerian cycle
    cycle_finder = CycleFinder(g, k, d)
    cycle = cycle_finder.find_cycle()

    # remove extra edge
    path = cycle

    # return path
    return path
