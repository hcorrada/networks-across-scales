import sys

from graph import build_graph
from path_finder import find_eulerian_path, find_valid_eulerian_path
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


# main assembler routine for single paths
# prints assembled string, or returns None if
# single path does not generate a valid string
def single_path(filename):
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

# create a function that checks path validity
def make_validity_function(k, d):
    # try to build string from path
    # return True if successful
    def check_path(path):
        string = build_string(path, k, d)
        return string is not None
    return check_path

# main assembler routine where all possible paths
# may be generated, it generates paths until a valid is one
# is found, it returns None if no valid path is found
def all_paths(filename):
    # read input from file
    paired_kmers, k, d = readdat(filename)

    # build DeBruijn graph from paired kmers
    g = build_graph(paired_kmers)

    # find Eulerian path that generates valid string
    # if there is no valid Eulerian path, this will return None
    path = find_valid_eulerian_path(g, make_validity_function(k,d))

    # check if a valid path was found
    if path is None:
        return None

    # get string from valid path
    string = build_string(path, k, d)
    print string

# this is here to play nicely with ipython
if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    all_paths(filename)
