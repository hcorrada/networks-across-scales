from graph import build_graph

# read input from file
# returns
#   paired_kmers: list of strings, each representing a pair of k-mers
#   k: k-mer length in pair
#   d: gap between k-mers in pair
def readdat(filename):
    with open(filename, 'r') as f:
        k, d = [int(s) for s in f.readline().strip().split()]
        paired_kmers = []
        for line in f:
            paired_kmers.append(line.strip())
        return paired_kmers, k, d

# read input from file
paired_kmers, k, d = readdat('test.txt')

# build DeBruijn graph from paired kmers
g = build_graph(paired_kmers)
print g.debug_print()
print
print g
