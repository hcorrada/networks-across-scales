import sys

def pattern_match(pattern, genome):
    k = len(pattern)
    ind = lambda x: slice(x, x+k)
    res = list()

    for i in xrange(len(genome) - k +1):
        if genome[ind(i)] == pattern:
            res.append(i)
    return res

filename = sys.argv[1]
with open(filename,'r') as f:
    pattern = f.readline().strip()
    genome = f.readline().strip()
    res = pattern_match(pattern, genome)
    print ' '.join(map(str, res))
