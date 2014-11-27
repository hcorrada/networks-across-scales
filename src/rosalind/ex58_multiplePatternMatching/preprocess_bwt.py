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


