import sys
import graph
from pathfinder import find_eulerian_path

# read input from file
def readdat(filename):
    with open(filename, 'r') as f:
        k, d = [int(s) for s in f.readline().strip().split()]
        paired_kmers = []
        for line in f:
            paired_kmers.append(line.strip())
        return paired_kmers, k, d

def build_graph(paired_kmers):
    return None

def main(filename):
    paired_kmers, k, d = readdat(filename)
    g = graph.build_graph(paired_kmers)
    print g

    path = find_eulerian_path(g, k, d)
    print path
    
    if path is None:
        return None

    string = path.get_string()
    print string

if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
