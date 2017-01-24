import sys

def hamming_distance(pattern, target):
    d = 0
    k = len(pattern)
    for i in xrange(k):
        if pattern[i] != target[i]:
            d += 1
    return d
    
# compute distance between pattern and string
# equal to minimum hamming distance between pattern and
# k-mers in string
def pattern_string_distance(pattern, string):
    n = len(string)
    k = len(pattern)
    min_dist = k+1

    for i in xrange(n-k+1):
        cur_dist = hamming_distance(pattern, string[slice(i, i+k)])
        if cur_dist < min_dist:
            min_dist = cur_dist
    return min_dist

filename = sys.argv[1]
with open(filename, 'r') as f:
    pattern = f.readline().strip()
    dna = [x.strip() for x in f.readline().strip().split()]
    print sum(map(lambda s: pattern_string_distance(pattern, s), dna))
