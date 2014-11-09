import sys
from collections import deque, defaultdict
import numpy as np
import heapq
from Bio.SubsMat.MatrixInfo import blosum62 as bl62

# wrapper to get score from
# hash representing upper triangular matrix
def get_score(x, y, mat):
  return mat[(x,y)] if (x,y) in mat else mat[(y,x)]

# class implementing a graph
# nodes are stored in hash, indexed by id
# each entry in node hash is of class Node (defined below)
# edges are stored in hash, indexed by (source_id, target_id)
# each entry edge hash is of class Edge (defined below)
class Graph:
    def __init__(self):
        # in-edges adjacency list
        self._inedges = defaultdict(list)
        self._outedges = defaultdict(list)

        # weights hash
        self._weights = {}
        
    # prints out each edge in the graph
    def __repr__(self):
        out = ""
        for target, sources in self._inedges.iteritems():
            for souce in sources:
                out += "%s -> %s\n" % (str(target), str(source))
        return out

    # add a weighted edge between two nodes with given ids
    # adds nodes if not in the graph already
    # see Node class below for edge representation
    def add_edge(self, source_id, target_id, weight):
        self._inedges[target_id].append(source_id)
        self._outedges[source_id].append(target_id)
        self._weights[(source_id, target_id)] = weight

    # find nodes with no incoming edges
    # in the graph
    def find_sources(self):
        sources = []
        for node, edges in self._inedges.iteritems():
            if len(edges) == 0:
                sources.append(node)
        return sources
        

    # find a topological ordering of the graph
    # returns a list of node ids
    def topological_ordering(self):
        # use this hash to mark edges that have been visited
        edge_marker = dict.fromkeys(self._weights.keys(), False)

        # result goes here
        order = []

        # start with nodes with no incoming edges
        # as candidates
        candidates = deque(self.find_sources())

        # while there are any candidates left
        while len(candidates) > 0:
            # add the next candidate to the topological ordering
            a = candidates.pop()
            order.append(a._id)

            # mark all the outgoing edges of the candidate
            for edge in a._outedges:
                edge_marker[edge.id()] = True

                # check if all incoming edges for the target have been visited
                if all([edge_marker[x.id()] for x in edge._target._inedges]):
                    # add target to candidates if so
                    candidates.append(edge._target)
        # check if any edge is left unmarked
        if any([not x for x in edge_marker.values()]):
            return None
        return order

    # compute the longest path between source and sink in the graph using
    # dynamic programming
    def longestpath(self, source, sink):
        # first, determine a topological order
        order = self.topological_ordering()
        if order is None:
            return 0, None

        print 'order done'
        # find where the source and sink are in this order
        source_index = order.index(source)
        sink_index = order.index(sink)

        # restrict topological order to source and sink
        order = order[slice(source_index, sink_index + 1)]
        n = len(order)

        # initialize length and backtrack table
        s = -float('inf') * np.ones(n)
        backtrack = np.zeros(n, np.dtype('int'))
        s[0] = 0.
        backtrack[0] = -1

        # now fill in table in order
        for i in xrange(1,n):
            print '.',
            # grab the next node in order
            node = self[order[i]]

            # see which of the incoming edges gives
            # the longest path from source to this node
            # use a heap to select longest incoming edge
            choices = []
            for edge in node._inedges:
                if edge._source._id not in order:
                    continue
                    
                j = order.index(edge._source._id)
                length = s[j] + edge._weight
                heapq.heappush(choices, (-length, j))
            if len(choices) > 0:
                # fill in with max choice
                s[i] = -choices[0][0]
                backtrack[i] = choices[0][1]

        # make the path from the backtrack table
        path = outputpath(backtrack, order, len(order) - 1, [])
        return int(s[-1]), path

# reconstruct the longest path from the backtrack table
def outputpath(backtrack, order, i, path):
    # use this iterative version to avoid
    # deep recursion (txs Guido!)
    while True:
        if i == 0:
            return [order[0]] + path
        i, path = backtrack[i], [order[i]] + path

# weighted edge structure
# _source and _target are of class Node        
class Edge:
    def __init__(self, source, target, weight):
        self._source = source
        self._target = target
        self._weight = weight

    def __repr__(self):
        return "%s:%d" % (self._target._id, self._weight)

    # return the id tuples used as id in the Graph edges hash
    def id(self):
        return (self._source._id, self._target._id)
    
# class representing a graph node            
class Node:
    def __init__(self, id):
        self._id = id
        self._outedges = []
        self._inedges = []

    # print node information
    def __repr__(self):
        children_strings = [str(edge) for edge in self._outedges]
        return "%s -> %s" % (self._id, ",".join(children_strings))

    # add an edge between this node
    # and given successor
    def add_outgoing_edge(self, edge):
        self._outedges.append(edge)

    # add an edge between this node
    # and given predecessor
    def add_incoming_edge(self, edge):
        self._inedges.append(edge)

    # return the indegree of the node
    def indegree(self):
        return len(self._inedges)

    # return the outdegree of the node
    def outdegree(self):
        return len(self._outedges)
    
def readdat(filename):
    with open(filename, 'r') as f:
        string1 = f.readline().strip()
        string2 = f.readline().strip()
    return string1, string2

# make the graph corresponding to the
# global alignment problem
def make_global_align_graph(v, w, sigma):
    nv = len(v)
    nw = len(w)

    graph = Graph()
    for i in xrange(nv+1):
        for j in xrange(nw+1):
            # add the '->' edge (gap in v)
            if i > 0:
                graph.add_edge((i-1,j), (i,j), -sigma)

            # add the diagonal edge (match v[i] and w[j])
            if i > 0 and j > 0:
                weight = get_score(v[i-1], w[j-1], bl62)
                graph.add_edge((i-1,j-1), (i,j), weight)

            # add the down edge (gap in w)
            if j > 0:
                graph.add_edge((i,j-1), (i,j), -sigma)
    return graph

def make_alignment(path, v, w):
    x,y = path[0]
    string1 = ''
    string2 = ''
    for i in xrange(1,len(path)):
        u,z = path[i]
        string1 += v[u - 1] if x != u else '-'
        string2 += w[z - 1] if y != z else '-'
        x,y = u,z
    return string1, string2

def main(filename):
    string1, string2 = readdat(filename)
    print len(string1)
    print len(string2)
    graph = make_global_align_graph(string1, string2, 5)
    print 'graph complete'
    return
    score, path = graph.longestpath((0,0), (len(string1), len(string2)))

    print 'got path'
    string1, string2 = make_alignment(path, string1, string2)    
    print score
    print string1
    print string2
    
if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
