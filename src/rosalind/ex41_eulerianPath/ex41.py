import sys
from collections import defaultdict

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
    def add_edge(self, source_id, target_id):
        if source_id not in self:
            self.add_node(source_id)
        if target_id not in self:
            self.add_node(target_id)
        # see Node class below for edge representation
        self[source_id].add_child(self[target_id])

    # add an edge from string representation
    # used in rosalind
    def add_edges_from_string(self, string):
        source, targets = string.split(" -> ")
        for target in targets.split(","):
            self.add_edge(source, target)

    # find a node in a given path
    # with an unused outgoing edge
    # path is assumed to be a list of node ids
    def find_start(self, path):
        # if path is empty, return arbitrary node
        if len(path) == 0:
            return self._nodes.values()[0]._id

        for nodeid in path:
            node = self._nodes[nodeid]
            # check if node has unused outgoing edge
            if node._num_used < len(node._children):
                # return id if that's the case
                return node._id
        # no possible start node found
        return None

    # walk given path (assumed to be a cycle),
    # starting at given node
    # path is a list of node ids
    def walk_path(self, startid, path):
        if len(path) == 0:
            return path

        # all we have to do is reorder the path list
        index = path.index(startid)
        return path[index:] + path[1:(index+1)]

    # find a new cycle starting at given node
    def extend_path(self, startid):
        path = [startid]
        curnode = self._nodes[startid]

        # check if current node has an unused outgoing edge
        while curnode._num_used < len(curnode._children):
            # it does, so use the next unused outgoing edge
            nextnode = curnode._children[curnode._num_used]
            # keep track of number of outgoing edges used
            curnode._num_used += 1

            # add the next node to the path
            path.append(nextnode._id)

            # keep walking
            curnode = nextnode
        return path

    # find an eulerian cycle 
    def eulerian_cycle(self):
        path = []

        # find a node to start the cycle
        startid = self.find_start(path)

        # if None, there is no available place to start
        while startid is not None:
            # walk current path (cycle) from this startid
            path = self.walk_path(startid, path)

            # extend it with a new cycle
            path = path[:-1] + self.extend_path(startid)

            # find the next place to start
            startid = self.find_start(path)
        return path

    # balance the graph
    def balance(self):
        # calculate out degree
        nedges_out = defaultdict(int)
        for id, node in self._nodes.iteritems():
            nedges_out[id] = len(node._children)
        nedges_out = sorted(nedges_out.items())

        # calculate in degree
        nedges_in = defaultdict(int)
        for id, node in self._nodes.iteritems():
            if id not in nedges_in:
                nedges_in[id] = 0
                
            for child in node._children:
                nedges_in[child._id] += 1
        nedges_in = sorted(nedges_in.items())

        # find source and sink
        source = None
        sink = None

        for i in xrange(len(nedges_out)):
            x, nout = nedges_out[i]
            y, nin = nedges_in[i]
            assert x == y

            if nout == nin:
                # this node is balanced                
                continue

            if nout < nin:
                # outdegree smaller than indegree                
                sink = x
            else:
                # outdegree smaller than indegree                                
                source = x

        # balance the graph adding edge between
        # sink and source
        self.add_edge(sink,source)

        # return the sink so we know how
        # walk the graph
        return sink

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
        children_ids = [str(child._id) for child in self._children]
        return "%s -> %s: %d" % (self._id, ",".join(children_ids), self._num_used)

    # add an edge between this node
    # and given child
    def add_child(self, child):
        self._children.append(child)
        
# read data from file
def readdat(filename):
    # initialize the graph
    graph = Graph()
    
    with open(filename, 'r') as f:
        for line in f:
            # add an edge from each line
            graph.add_edges_from_string(line.strip())
    return graph


def main(filename):
    # read input and construct graph
    graph = readdat(filename)

    # balance the graph
    # by adding an edge
    sink = graph.balance()

    # find the eulerian cycle
    cycle = graph.eulerian_cycle()

    # remove the balancing edge
    # to turn cycle into a path
    index = cycle.index(sink)
    path = cycle[index+1:] + cycle[1:index+1]
    print "->".join(map(str,path))
    
if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
