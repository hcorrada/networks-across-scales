import sys

# compute the kmer composition of string
# input:
#   k: k-mer length
#   text: dna string
# output:
#   list of k-mers in text
def kmer_composition(k, text):
    kmers = []
    for i in xrange(len(text) - k +1):
        kmers.append(text[slice(i,i+k)])
    kmers.sort()
    return kmers

# read input from file
# expected format:
#   k
#   text
# output:
#   tuple k, text
def readdat(filename):
    with open(filename, 'r') as f:
        k = int(f.readline().strip())
        test = f.readline().strip()
    return k, test

# make debruijn graph, i.e., k-mer overlap graph
# input:
#   kmers: list of kmers
# output:
#   dictionary:
#       keys: k-1 mers (edge source)
#       values: k-mers (edge targets)
def makeGraph(kmers):
    graph = {}
    nkmers = len(kmers)
    for i in xrange(nkmers):
        kmer = kmers[i]
        source = kmer[:-1]
        target = kmer[1:]

        # if source is already a key in dictionary
        if source in graph:
            # append target to target list
            graph[source].append(target)
        else:
            # add to graph dictionary
            graph[source] = [target]
    # sort by source and return
    return sorted(graph.items())

def main(filename):
    k, text = readdat(filename)
    kmers = kmer_composition(k, text)
    graph = makeGraph(kmers)

    # print out adjacency list for each node
    out = [x + " -> " + ",".join(y) for (x,y) in graph]
    print "\n".join(out)

if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
