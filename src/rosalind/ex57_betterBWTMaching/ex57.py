import sys
import numpy as np
from collections import defaultdict

# tag each symbol with its occurences
# charlist is a list of characters
# returns a list of tuples, each tuple
# contains the character in the list, and its
# number occurence in the list
def count_occurences(charlist):
    counter = defaultdict(int)
    for c in charlist:
        counter[c] += 1
    symbols = counter.keys()

    count = np.zeros((len(charlist)+1, len(symbols)), dtype=np.int32)

    for i in xrange(len(charlist)):
        symbol = charlist[i]
        vec = np.array(count[i,])
        index = symbols.index(symbol)
        vec[index] += 1
        count[i+1,:] = vec

    def fn(symbol, position):
        index = symbols.index(symbol)
        return count[position, index]
    return fn

# prepocesses the bwt by tagging occurences in
# the bwt (the last column) and its sorted
# version (the first column)
# bwt is a string containing the bwt of some text
# returns first and last columns tagged with occurecences
def preprocess_bwt(bwt):
    bwt = list(bwt)
    count = count_occurences(bwt)

    bwt = sorted(bwt)
    first_occurence = {}
    for i in xrange(len(bwt)):
        c = bwt[i]
        if c not in first_occurence:
            first_occurence[c] = i 

    return first_occurence, count

# counts the number of times pattern occurs in text using
# first and last columns of M(text) matrix reconstructed
# using preprocess_bwt
def count_matches(pattern, first_occurence, count, n):
    # have top and bottom pointers at first and last
    # row of M matrix
    top = 0
    bottom = n - 1

    # while there are rows to be checked
    while top <= bottom:
        # and patterns symbols to be checked
        if len(pattern) > 0: 
            symbol = pattern[-1] # check the last symbol in pattern
            pattern = pattern[:-1] # remove last symbol in pattern

            top = first_occurence[symbol] + count(symbol, top)
            bottom = first_occurence[symbol] + count(symbol, bottom + 1) - 1
        else:
            return bottom - top + 1
    return 0

def readdat(filename):
    with open(filename, 'r') as f:
        bwt = f.readline().strip()
        patterns = [x.strip() for x in f.readline().split()]
        return bwt, patterns

def main(filename):
    bwt, patterns = readdat(filename)
    first_occurence, count = preprocess_bwt(bwt)

    out = []
    for pattern in patterns:
        out.append(count_matches(pattern, first_occurence, count, len(bwt)))
    print ' '.join(map(str,out))

if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
