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
        self._available_nodes = set()

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
            next_node = self.next_target(cur_node)
            if self.is_available(cur_node):
                self._available_nodes.add(cur_node)
            elif cur_node in self._available_nodes:
                self._available_nodes.discard(cur_node)
            cur_node = next_node
        path.append(cur_node)
        return path

    def run(self):
        cycle = Path()
        node = iter(self._graph).next()
        self._available_nodes.add(node)

        while len(self._available_nodes) > 0:
            node = self._available_nodes.pop()
            #print "next node"
            #print node

            item = cycle.find_item(node)
            new_cycle = self.make_a_cycle(node)
            #print new_cycle

            cycle = cycle.stitch(item, new_cycle)
            #print cycle
            #print
        return cycle
