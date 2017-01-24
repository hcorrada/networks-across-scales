import sys

def overlapGraph(kmers):
    k = len(kmers[0])
    graph = {}
    for i in xrange(len(kmers)):
        kmer = kmers[i]
        graph[kmer] = []
        for j in xrange(len(kmers)):
            if j == i:
                continue
            if kmers[j][:-1] == kmer[1:]:
                graph[kmer].append(kmers[j])

    adjList = []
    for kmer, children in graph.iteritems():
        for kmer2 in children:
            adjList.append((kmer,kmer2))
    adjList.sort()
    return adjList

def readdat(filename):
    with open(filename, 'r') as f:
        kmers = [line.strip() for line in f]
        return kmers
    
def main(filename):
    kmers = readdat(filename)
    graph = overlapGraph(kmers)
    print "\n".join([x + " -> " + y for (x,y) in graph])
    
if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
