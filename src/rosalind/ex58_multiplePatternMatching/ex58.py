import sys
from construct_bwt import BWT
from suffix_array import SuffixArray, PartialSuffixArray

def readdat(filename):
    with open(filename, 'r') as f:
        text = f.readline().strip()
        patterns = [x.strip() for x in f.readlines()]
        return text, patterns


def main(filename):
    k = 5
    
    text, patterns = readdat(filename)
    text = text + '$'
    
    bwt = BWT(text)
    sa = PartialSuffixArray(text, k, bwt)
    
    indices = []
    for pattern in patterns:
        top, bottom = bwt.find_matches(pattern)
        if top != -1:
            indices += sa.get_match_indices(top, bottom)
    print ' '.join(map(str, sorted(indices)))

if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
