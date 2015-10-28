import sys

from graph import build_graph
from path_finder import find_eulerian_path
from string_builder import build_string

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


# main assembler routine
# prints assembled string
def main(filename):
    # read input from file
    paired_kmers, k, d = readdat(filename)

    # build DeBruijn graph from paired kmers
    g = build_graph(paired_kmers)

    # find Eulerian path in Debruijn graph
    # if there is no Eulerian path, this will return None
    path = find_eulerian_path(g)

    # check if an Eulerian path was found
    if path is None:
        return None

    # if so, construct a string from the given path
    # if path does not spell out a valid string
    # this will return None
    string = build_string(path, k, d)

    # print resulting string
    print string

# this is here to play nicely with ipython
if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
