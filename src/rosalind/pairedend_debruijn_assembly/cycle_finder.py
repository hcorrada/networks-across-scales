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
    def __init__(self, g):
        self._graph = g

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
        path = Path()

        cur_node = start_node
        while self.is_available(cur_node):
            path.append(cur_node)
            cur_node = self.next_target(cur_node)
        return path

    def run(self):
        cycle = Path()

        for node in self._graph:
            if not self.is_available(node):
                continue

            print "next node:"
            print node

            item = cycle.find_item(node)
            print item
            
            new_cycle = self.make_a_cycle(node)
            print new_cycle
            cycle = cycle.stitch(item, new_cycle)
            print cycle
            print

        return cycle
