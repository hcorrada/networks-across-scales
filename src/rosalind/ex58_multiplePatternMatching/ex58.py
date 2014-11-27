import sys
from construct_bwt import get_bwt
from bwt_match import preprocess_bwt, find_matches
from suffix_array import make_suffix_array

def readdat(filename):
    with open(filename, 'r') as f:
        text = f.readline().strip()
        patterns = [x.strip() for x in f.readlines()]
        return text, patterns

def main(filename):
    text, patterns = readdat(filename)
    text = text + '$'
    bwt = get_bwt(text)
    first_occurence, count = preprocess_bwt(bwt)

    sa = make_suffix_array(text)
    
    indices = []
    for pattern in patterns:
        top, bottom = find_matches(pattern, first_occurence, count, len(bwt))
        if top != -1:
            for i in xrange(top, bottom+1):
                indices.append(sa[i])
    print ' '.join(map(str, sorted(indices)))

if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
