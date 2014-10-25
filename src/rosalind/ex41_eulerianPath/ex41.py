import sys
from collections import defaultdict

class Graph:
    def __init__(self):
        self._nodes = {}

    def __repr__(self):
        out = ""
        for node in self._nodes.itervalues():
            out += node.__repr__() + "\n"
        return out
        
    def __getitem__(self, i):
        return self._nodes[i]

    def __contains__(self, i):
        return i in self._nodes

    def keys(self):
        return self._nodes.keys()
        
    def add_node(self, id):
        self._nodes[id] = Node(id)

    def add_edge(self, source_id, target_id):
        if source_id not in self:
            self.add_node(source_id)
        if target_id not in self:
            self.add_node(target_id)
        self[source_id].add_child(self[target_id])

    def add_edges_from_string(self, string):
        source, targets = string.split(" -> ")
        for target in targets.split(","):
            self.add_edge(source, target)

    def find_start(self, path):
        if len(path) == 0:
            return self._nodes.values()[0]._id

        for nodeid in path:
            node = self._nodes[nodeid]
            if node._num_used < len(node._children):
                return node._id
        return None

    def walk_path(self, startid, path):
        if len(path) == 0:
            return path
        
        index = path.index(startid)
        return path[index:] + path[1:(index+1)]

    def extend_path(self, startid):
        path = [startid]
        curnode = self._nodes[startid]
        while curnode._num_used < len(curnode._children):
            nextnode = curnode._children[curnode._num_used]
            curnode._num_used += 1
            path.append(nextnode._id)
            curnode = nextnode
        return path
    
    def eulerian_cycle(self):
        path = []
        startid = self.find_start(path)
        while startid is not None:
            path = self.walk_path(startid, path)
            path = path[:-1] + self.extend_path(startid)
            startid = self.find_start(path)
        return path

    def balance(self):
        nedges_out = defaultdict(int)
        for id, node in self._nodes.iteritems():
            nedges_out[id] = len(node._children)
        nedges_out = sorted(nedges_out.items())

        nedges_in = defaultdict(int)
        for id, node in self._nodes.iteritems():
            if id not in nedges_in:
                nedges_in[id] = 0
                
            for child in node._children:
                nedges_in[child._id] += 1
        nedges_in = sorted(nedges_in.items())

        source = None
        sink = None
        
        for i in xrange(len(nedges_out)):
            x, nout = nedges_out[i]
            y, nin = nedges_in[i]
            assert x == y
            if nout == nin:
                continue
            if nout < nin:
                sink = x
            else:
                source = x
        self.add_edge(sink,source)
        return sink
        
class Node:
    def __init__(self, id):
        self._id = id
        self._children = []
        self._num_used = 0

    def __repr__(self):
        children_ids = [str(child._id) for child in self._children]
        return "%s -> %s: %d" % (self._id, ",".join(children_ids), self._num_used)
    
    def add_child(self, child):
        self._children.append(child)
        

def readdat(filename):
    graph = Graph()
    
    with open(filename, 'r') as f:
        for line in f:
            graph.add_edges_from_string(line.strip())
    return graph

def main(filename):
    graph = readdat(filename)
    sink = graph.balance()
    cycle = graph.eulerian_cycle()
    index = cycle.index(sink)
    path = cycle[index+1:] + cycle[1:index+1]
    print "->".join(map(str,path))
    
if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
