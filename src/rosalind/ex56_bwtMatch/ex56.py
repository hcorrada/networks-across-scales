import sys

import sys
from collections import defaultdict

# tag each symbol with its occurences
# charlist is a list of characters
# returns a list of tuples, each tuple
# contains the character in the list, and its
# number occurence in the list
def tag_occurences(charlist):
    counter = defaultdict(int)
    out = []
    for c in charlist:
        counter[c] += 1
        out.append((c,counter[c]))
    return out

# prepocesses the bwt by tagging occurences in
# the bwt (the last column) and its sorted
# version (the first column)
# bwt is a string containing the bwt of some text
# returns first and last columns tagged with occurecences
def preprocess_bwt(bwt):
    bwt = list(bwt)
    last_column = tag_occurences(bwt)    
    first_column = tag_occurences(sorted(bwt))
    return first_column, last_column

# counts the number of times pattern occurs in text using
# first and last columns of M(text) matrix reconstructed
# using preprocess_bwt
def count_matches(pattern, first_column, last_column):
    # have top and bottom pointers at first and last
    # row of M matrix
    top = 0
    bottom = len(first_column)-1

    # while there are rows to be checked
    while top <= bottom:
        # and patterns symbols to be checked
        if len(pattern) > 0: 
            symbol = pattern[-1] # check the last symbol in pattern
            pattern = pattern[:-1] # remove last symbol in pattern

            # find the indices of first and last occurences of symbol in current range of rows
            # in the last column
            topIndex,  bottomIndex = find_symbol(symbol, last_column[top:(bottom+1)])
            if topIndex == -1:
                # symbol not found, so pattern does not occur in text
                return 0

            # now find the new range of rows to be checked with next symbol
            bottom = first_to_last(top+ bottomIndex, first_column, last_column)
            top = first_to_last(top + topIndex, first_column, last_column)
            
        else:
            return bottom - top + 1
    return 0

# find the row in first column corresponding to character
# in last column at row index
def first_to_last(index, first_column, last_column):
    return first_column.index(last_column[index])

# find symbol in last column
# returns indices of first and last occurence of symbol in the column
def find_symbol(symbol, last_column):
    # extract the symbols 
    symbols, _ignored = zip(*last_column)

    # symbol not there, so return -1
    if symbol not in symbols:
        return -1, -1

    # find the first occurence
    firstIndex = symbols.index(symbol)

    # reverse the symbols
    symbols = list(reversed(symbols))

    # find the last occurence (first occurence in reversed string)
    lastIndex = len(symbols) - symbols.index(symbol) - 1
    return firstIndex, lastIndex

def readdat(filename):
    with open(filename, 'r') as f:
        bwt = f.readline().strip()
        patterns = [x.strip() for x in f.readline().split()]
        return bwt, patterns

def main(filename):
    bwt, patterns = readdat(filename)
    first_column, last_column = preprocess_bwt(bwt)

    out = []
    for pattern in patterns:
        out.append(count_matches(pattern, first_column, last_column))
    print ' '.join(map(str,out))

if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
