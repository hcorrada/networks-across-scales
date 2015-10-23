import sys

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
    graph = build_graph(paired_kmers)
    path = find_eulerian_path(paired_kmers, k, d)
    string = path.get_string()
    print string

if __name__ == '__main__' and and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
