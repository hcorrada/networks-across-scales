import sys

def kmer_composition(k, text):
    kmers = []
    for i in xrange(len(text) - k +1):
        kmers.append(text[slice(i,i+k)])
    kmers.sort()
    return kmers

def readdat(filename):
    with open(filename, 'r') as f:
        k = int(f.readline().strip())
        test = f.readline().strip()
    return k, test

def makeGraph(k, kmers):
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
    
def main(filename):
    k, text = readdat(filename)
    kmers = kmer_composition(k, text)
    graph = makeGraph(k, kmers)
    out = [x + " -> " + ",".join(y) for (x,y) in graph]
    print "\n".join(out)

if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
