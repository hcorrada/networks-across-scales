import sys
from bwt import BWT
from suffix_array import SuffixArray, PartialSuffixArray
from itertools import izip_longest

def readdat(filename):
    with open(filename, 'r') as f:
        text = f.readline().strip()
        patterns = [x.strip() for x in f.readline().split()]
        d = int(f.readline().strip())
        return text, patterns, d

def make_seeds(pattern, d):
    n = len(pattern)
    k = n / (d + 1)
    args = [iter(pattern)] * k
    seeds = [''.join(x) for x in izip_longest(fillvalue='', *args)]
    if len(seeds) > d+1:
        last = seeds[-1]
        seeds = seeds[:-1]
        seeds[-1] += last
    return seeds, k

def find_approximate_candidates(pattern, bwt, sa, d):
    seeds, k = make_seeds(pattern, d)
    indices = []
    for i in xrange(len(seeds)):
        seed = seeds[i]
        top, bottom = bwt.find_matches(seed)
        if top == -1:
            continue
        candidates = sa.get_match_indices(top, bottom)
        for index in candidates:
            candidate = index - i*k
            if candidate not in indices:
                indices.append(index - i*k)
    return sorted(indices)

def within_edit_distance(x, y, d):
    dist = 0
    for i in xrange(len(x)):
        dist += 1 if x[i] != y[i] else 0
        if dist > d:
            return False
    return True

def find_approximate_matches(pattern, bwt, sa, d, text):
    candidates = find_approximate_candidates(pattern, bwt, sa, d)
    candidates = filter(lambda x: within_edit_distance(pattern, text[x:], d), candidates)
    return candidates

def main(filename):
    k = 5
    c = 5
    
    text, patterns, d = readdat(filename)
    text = text + '$'

    bwt = BWT(text, checkpoints=c)
    sa = PartialSuffixArray(text, k, bwt)
    
    indices = []
    for pattern in patterns:
        indices += find_approximate_matches(pattern, bwt, sa, d, text)
    print ' '.join(map(str, sorted(indices)))

if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
