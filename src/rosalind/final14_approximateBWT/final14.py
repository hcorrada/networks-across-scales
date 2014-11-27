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

def main(filename):
    k = 5
    c = 5
    
    text, patterns, d = readdat(filename)
    text = text + '$'

    bwt = BWT(text, checkpoints=c)
    sa = PartialSuffixArray(text, k, bwt)
    
    indices = []
    for pattern in patterns[:1]:
        top, bottom = bwt.find_matches(pattern)
        if top != -1:
            indices += sa.get_match_indices(top, bottom)
    print ' '.join(map(str, sorted(indices)))

if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
