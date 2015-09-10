import sys

def find_pattern(pattern, genome):
    n = len(genome)
    k = len(pattern)
    res = list()

    for i in xrange(n - k + 1):
        if genome[slice(i, i + k)] == pattern:
            res.append(i)
    return res

filename = sys.argv[1]
with open(filename,'r') as f:
    pattern = f.readline().strip()
    genome = f.readline().strip()
    res = find_pattern(pattern, genome)
    print ' '.join(map(str, res))
