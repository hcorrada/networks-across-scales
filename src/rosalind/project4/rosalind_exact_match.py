import sys
from approximate_matcher.bwt import BWT, _preprocess_bwt

def readdat(filename):
    with open(filename, 'r') as f:
        bwt = f.readline().strip()
        patterns = [x.strip() for x in f.readline().split()]
        return bwt, patterns

# read target, patterns and number of mismatches
def main(filename):
    bwt, patterns = readdat(filename)

    # initialize approximate matcher
    bwt_obj = BWT()
    bwt_obj._bwt = bwt
    bwt_obj.first_occurence, bwt_obj.count = _preprocess_bwt(bwt)

    # find approximate matches for each pattern
    result = []
    for pattern in patterns:
        top, bottom = bwt_obj._get_matching_rows(pattern)
        result.append(0 if top == -1 else bottom - top + 1)
    print ' '.join(map(str, result))

# play nicely with ipython
if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
