import sys

from graph import build_graph
from path_finder import find_eulerian_path
from string_builder import build_string

# read input from file
def readdat(filename):
    with open(filename, 'r') as f:
        k, d = [int(s) for s in f.readline().strip().split()]
        paired_kmers = []
        for line in f:
            paired_kmers.append(line.strip())
        return paired_kmers, k, d


def main(filename):
    paired_kmers, k, d = readdat(filename)
    g = build_graph(paired_kmers)
    print g

    path = find_eulerian_path(g)
    print path

    if path is None:
        return None

    string = build_string(path, k, d)
    print string

if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
