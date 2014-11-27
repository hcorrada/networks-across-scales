import numpy as np
from collections import defaultdict

# precompute the count table
# return a function that returns the
# precomputed count
def count_occurences(charlist):
    # figure out symbols occuring in text
    counter = defaultdict(int)
    for c in charlist:
        counter[c] += 1
    symbols = counter.keys()

    # allocate count table
    count = np.zeros((len(charlist)+1, len(symbols)), dtype=np.int32)

    # fill in count table with running count
    # of occurence counts
    for i in xrange(len(charlist)):
        symbol = charlist[i]
        vec = np.array(count[i,])
        index = symbols.index(symbol)
        # increase running count for symbol
        vec[index] += 1
        count[i+1,:] = vec

    # function that return precomputed count
    def fn(symbol, position):
        index = symbols.index(symbol)
        return count[position, index]
    return fn

# precompute the count table but only store
# rows that are multiple of 'checkpoints' argument
# returns a function that returns the count
def count_occurences_checkpoints(charlist, checkpoints):
    # figure out symbols occuring in text
    counter = defaultdict(int)
    for c in charlist:
        counter[c] += 1
    symbols = counter.keys()

    # allocate checkpoint table
    nrows = (len(charlist)-1) / checkpoints 
    count = np.zeros((nrows+1, len(symbols)), dtype=np.int32)

    # vector of running counts
    vec = np.array(count[0,:])

    # update running counts
    # and fill in checkpoint table 
    for i in xrange(len(charlist)):
        # this is a checkpoint row
        if (i % checkpoints) == 0:
            rowindex = i / checkpoints
            count[rowindex, :] = np.array(vec)
        # update running count
        symbol = charlist[i]
        index = symbols.index(symbol)
        vec[index] += 1
    # check the last entry
    if (i % checkpoints) == 0:
        rowindex = i / checkpoints
        count[rowindex,:] = np.array(vec)

    # function that returns count
    def fn(symbol, position):
        # find column in checkpoint table
        index = symbols.index(symbol)
        # find row in checkpoint table
        rowindex = position / checkpoints
        rowindex = min(rowindex, nrows)

        # get count from checkpoint row
        out = count[rowindex, index]

        # now update count by walking bwt
        # find checkpoint position in bwt
        index = rowindex * checkpoints

        # update count up to required position
        while index < position:
            # update count if matching symbol
            out += 1 if charlist[index] == symbol else 0
            index += 1
        return out

    # return the function
    return fn

# return the first occurence of each symbol
# in the first column of rotation matrix
def get_first_occurence(bwt):
    bwt = sorted(bwt) # sort bwt
    first_occurence = {}

    # find first occurences
    for i in xrange(len(bwt)):
        c = bwt[i]
        if c not in first_occurence:
            first_occurence[c] = i 
    return first_occurence



