import sys

def debruijnGraph(kmers):
    k = len(kmers[0])
    graph = {}
    for i in xrange(len(kmers)):
        kmer = kmers[i]
        source = kmer[:-1]
        target = kmer[1:]
        if source in graph:
            graph[source].append(target)
        else:
            graph[source] = [target]
    return sorted(graph.items())

def readdat(filename):
    with open(filename, 'r') as f:
        kmers = [line.strip() for line in f]
        return kmers

def main(filename):
    kmers = readdat(filename)
    graph = debruijnGraph(kmers)
    print "\n".join([x + " -> " + ",".join(y) for (x,y) in graph])
    
if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
