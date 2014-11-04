import sys
from collections import deque
import numpy as np
import heapq

# class implementing a graph
# nodes are stored in hash, indexed by id
# each entry in hash is of class Node (defined below)
class Graph:
    def __init__(self):
        self._nodes = {}

    # prints out each node in the graph
    def __repr__(self):
        out = ""
        for node in self._nodes.itervalues():
            out += node.__repr__() + "\n"
        return out

    # returns node by id, so graph[id] works
    def __getitem__(self, i):
        return self._nodes[i]

    # checks if a node is in the graph already 'id in graph'
    def __contains__(self, i):
        return i in self._nodes

    # returns the list of node ids
    def keys(self):
        return self._nodes.keys()

    # add a new node with given id
    def add_node(self, id):
        self._nodes[id] = Node(id)

    # adds an edge between two nodes with given ids
    # adds nodes if not in the graph already
    # see Node class below for edge representation
    def add_edge(self, source_id, target_id, weight):
        if source_id not in self:
            self.add_node(source_id)
        if target_id not in self:
            self.add_node(target_id)
        # see Node class below for edge representation
        self[source_id].add_child(self[target_id], weight)

    # add an edge from string representation
    # used in rosalind
    def add_edge_from_string(self, string):
        source, tmp = string.split("->")
        target, weight = tmp.split(":")
        self.add_edge(source, target, int(weight))

    def indegree(self, x):
        cnt = 0
        for node in self._nodes.itervalues():
            for edge in node._children[node._num_used:]:
                if edge._target == x:
                    cnt += 1
        return cnt

    def any_edges_unused(self):
        for node in self._nodes.itervalues():
            if node._num_used < len(node._children):
                return True
        return False

    def find_sources(self):
        nonsources = set()
        for node in self._nodes.itervalues():
            for edge in node._children:
                nonsources.add(edge._target._id)
        sources = []
        for node in self._nodes.itervalues():
            if node._id not in nonsources:
                sources.append(node)
        return sources
        
        
    def topological_ordering(self):
        order = []
        candidates = deque(self.find_sources())
        while len(candidates) > 0:
            a = candidates.pop()
            order.append(a._id)
            for edge in a._children:
                a._num_used += 1
                if self.indegree(edge._target) == 0:
                    candidates.append(edge._target)
        if self.any_edges_unused():
            return None
        return order
    
    def longestpath(self, source, sink):
        order = self.topological_ordering()
        if order is None:
            return 0, None
        
        source_index = order.index(source)
        sink_index = order.index(sink)
        order = order[slice(source_index, sink_index + 1)]
        n = len(order)

        s = -float('inf') * np.ones(n)
        backtrack = np.zeros(n, np.dtype('int'))
        s[0] = 0.
        backtrack[0] = -1
        
        for i in xrange(1,n):
            node = self[order[i]]
            choices = []
            for j in xrange(i):
                candidate = self[order[j]]
                edge = candidate.get_edge(node._id)
                if edge is not None:
                    length = s[j] + edge._weight
                    heapq.heappush(choices, (-length, j))
            if len(choices) > 0:
                s[i] = -choices[0][0]
                backtrack[i] = choices[0][1]

        path = outputpath(backtrack, order, len(order) - 1, [])
        return int(s[-1]), path

def outputpath(backtrack, order, i, path):
    while True:
        if i == 0:
            return [order[0]] + path
        i, path = backtrack[i], [order[i]] + path

class Edge:
    def __init__(self, target, weight):
        self._target = target
        self._weight = weight

    def __repr__(self):
        return "%s:%d" % (self._target._id, self._weight)
    
# class representing a graph node            
class Node:
    def __init__(self, id):
        self._id = id
        self._children = []
        # this keeps track of number
        # of outgoing edges used so far
        self._num_used = 0

    # print node information
    def __repr__(self):
        children_strings = [str(edge) for edge in self._children]
        return "%s -> %s" % (self._id, ",".join(children_strings))

    # add an edge between this node
    # and given child
    def add_child(self, child, weight):
        edge = Edge(child, weight)
        self._children.append(edge)

    def get_edge(self, id):
        for edge in self._children:
            if edge._target._id == id:
                return edge
        return None
            
def readdat(filename):
    with open(filename, 'r') as f:
        source = f.readline().strip()
        sink = f.readline().strip()

        graph = Graph()
        
        for line in f:
            graph.add_edge_from_string(line.strip())

    return source, sink, graph

def main(filename):
    source, sink, graph = readdat(filename)
    weight, path = graph.longestpath(source, sink)
    print weight
    print "->".join(path)
    
if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
