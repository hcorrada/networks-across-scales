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

    def find_source(self):
        if len(self._nodes) == 0:
            return None
        
        for node in self._nodes.values():
            if node.get_outdegree() > node.get_indegree():
                return node
        return None

    def reconstruct(self, path):
        edge = path[0]
        text = edge[0]
        for edge in path[1:]:
            kmer = edge[0]
            text += kmer[-1]
        kmer = edge[1]
        text += kmer[-1]
        return text
    
    def get_next_contig(self):
        source = self.find_source()
        if source is None:
            return None, 0
        target, numedges = source.pop_next_child(None)
        path = [(source._id, target._id)]
        numcopies = numedges
            
        while True:
            source = target
            if source._sink or source.get_numchildren() > 1:
                break
            
            target, numedges = source.peek_next_child()
            if numedges < numcopies or target.get_numchildren() > 1:
                break
            
            if numedges > numcopies:
                break
            
            target, _ignored = source.pop_next_child(numcopies)
            edge = (source._id, target._id)
            path.append(edge)

        contig = self.reconstruct(path)
        self.prune()
        return contig, numcopies
            
    def prune(self):
        pruned = []
        for node_id, node in self._nodes.iteritems():
            if node.get_indegree() == 0 and node.get_outdegree() == 0:
                pruned.append(node_id)
        for node_id in pruned:
            del self._nodes[node_id]
                
    def get_contigs(self):
        res = []
        contig, numcopies = self.get_next_contig()

        while contig is not None:
            res += [contig] * numcopies
            contig, numcopies = self.get_next_contig()
        return res

    def flag_sinks(self):
        for node in self._nodes.itervalues():
            if node.get_indegree() > node.get_outdegree():
                node._sink = True
    
class Node:
    def __init__(self, id):
        self._id = id
        self._children = defaultdict(int)
        self._outdegree = 0
        self._indegree = 0
        self._sink = False

    def __repr__(self):
        children_strings = ["%s (%d)" % (str(child._id), numedges) for child, numedges in self._children.iteritems()]
        return "%s -> %s: %d - %d : %s" % (self._id, ",".join(children_strings), self.get_indegree(), self.get_outdegree(), "sink" if self._sink else "")

    def get_numchildren(self):
        return sum([x > 0 for x in self._children.values()])
    
    def get_outdegree(self):
        return sum(self._children.values())

    def get_indegree(self):
        return self._indegree

    def add_child(self, child):
        self._children[child] += 1
        child._indegree += 1

    def peek_next_child(self):
        for child, numedges in self._children.iteritems():
            if numedges > 0:
                return child, numedges
        return None, 0
            
    def pop_next_child(self, numcopies):        
        child, numedges = self.peek_next_child()
        if child is None:
            return None
        
        if numcopies is None:
            numcopies = numedges
        self._children[child] -= numcopies
        child._indegree -= numcopies
        return child, numedges

def debruijnGraph(kmers):
    k = len(kmers[0])
    graph = Graph()
    
    for i in xrange(len(kmers)):
        kmer = kmers[i]
        source = kmer[:-1]
        target = kmer[1:]
        graph.add_edge(source, target)
    return graph

def readdat(filename):
    with open(filename, 'r') as f:
        kmers = [line.strip() for line in f]
        return kmers
                
def main(filename):
    kmers = readdat(filename)
    graph = debruijnGraph(kmers)
    graph.flag_sinks()
    print graph
    contigs = graph.get_contigs()
    print " ".join(contigs)
    
if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
