import sys

from approximate_matcher import ApproximateMatcher

def readdat(filename):
    with open(filename, 'r') as f:
        text = f.readline().strip()
        patterns = [x.strip() for x in f.readline().split()]
        d = int(f.readline().strip())
        return text, patterns, d

# read target, patterns and number of mismatches
def main(filename):
    text, patterns, d = readdat(filename)

    # initialize approximate matcher
    am = ApproximateMatcher(text)

    # find approximate matches for each pattern
    indices = []
    for pattern in patterns:
        indices += am.get_matches(pattern, d)
    print ' '.join(map(str, sorted(indices)))

# play nicely with ipython
if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
