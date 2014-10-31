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
    
    def eulerian_cycle(self, startid=None):
        path = []
        if startid is None:
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
        return source, sink

    def eulerian_path(self):
        source, sink = self.balance()
        cycle = self.eulerian_cycle(sink)
        index = cycle.index(sink)
        path = cycle[index+1:] + cycle[1:index+1]
        return path
        
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

def debruijnGraph(kmers):
    k = len(kmers[0])
    graph = Graph()
    
    for i in xrange(len(kmers)):
        kmer = kmers[i]
        read1, read2 = kmer.split("|")
        source1 = read1[:-1]
        source2 = read2[:-1]
        
        target1 = read1[1:]
        target2 = read2[1:]
        
        source = source1 + "|" + source2
        target = target1 + "|" + target2
        
        graph.add_edge(source, target)
    return graph

def readdat(filename):
    with open(filename, 'r') as f:
        d = int(f.readline().strip())
        kmers = [line.strip() for line in f]
        return d, kmers

def reconstruct(path, d):
    read_pair = path[0]
    text1, text2 = read_pair.split("|")
    k = len(text1) + 1
    
    for read_pair in path[1:]:
        kmer1, kmer2 = read_pair.split("|")
        text1 += kmer1[-1]
        text2 += kmer2[-1]
    return text1 + text2[-(k+d):]
    
def main(filename):
    d, kmers = readdat(filename)
    graph = debruijnGraph(kmers)
    path = graph.eulerian_path()
    text = reconstruct(path, d)
    print text
    
if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
