from path import Path
from collections import defaultdict
import itertools

class ProbedIter:
    def __init__(self, it):
        self._it, self._next_it = itertools.tee(iter(it))
        self._advance()

    def _advance(self):
        try:
            self._lookahead = self._next_it.next()
        except StopIteration:
            self._lookahead = None

    def next(self):
        self._advance()
        return self._it.next()

    def has_more(self):
        return self._lookahead is not None

class CycleFinder:
    def __init__(self, g, k, d):
        self._graph = g
        self._k = k
        self._d = d

        self._target_iterators = dict([(node.label(), ProbedIter(node.target_labels())) for node in self._graph])
        self._cycle_items = defaultdict(list)

    def next_target(self, node):
        it = self._target_iterators[node.label()]
        if it.has_more():
            label = it.next()
            return self._graph[label]
        else:
            return None

    def is_available(self, node):
        it = self._target_iterators[node.label()]
        return it.has_more()

    def make_a_cycle(self, start_node):
        path = Path(self._k, self._d)

        cur_node = start_node
        while self.is_available(cur_node):
            path.append(cur_node)
            cur_node = self.next_target(cur_node)
        return path

    def find_cycle(self):
        cycle = Path(self._k, self._d)

        for node in self._graph:
            if not self.is_available(node):
                continue

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
