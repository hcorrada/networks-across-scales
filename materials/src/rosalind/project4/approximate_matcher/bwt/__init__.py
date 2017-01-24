from preprocess_bwt import _get_first_occurence_fn, _get_count_fn

# THIS IS A STUB, YOU NEED TO IMPLEMENT THIS
#
# Construct the Burrows-Wheeler transform for given text
# also compute the suffix array
#
# Input:
#   text: a string (character `$` assumed to be last character)
#
# Output:
#   a tuple (bwt, suffix_array):
#       bwt: string containing the Burrows-Wheeler transform of text
#       suffix_array: the suffix array of text
def _construct(text):
    # rotate text so that last character of rotation is on the appropriate position
    # i.e. the ith position has the last character of the rotation that
    # starts at the ith position in text (this is not the same order that
    # the rotations are listed on pg. 326 of Compeau & Pevzner
    rotated_text = text[-1] + text[:-1]

    # make tuples of character in bwt and the rotation it corresponds to
    indexed_chars = zip(list(rotated_text), [i for i in xrange(len(text))])

    # sort tuples by corresponding rotations
    # the key function return the rotated text starting at appropriate index
    sorted_chars = sorted(indexed_chars, key=lambda x: text[x[1]:] + text[:x[1]])

    # extract the bwt characters after sorting rotations
    bwt_chars, indices = zip(*sorted_chars)
    bwt = ''.join(list(bwt_chars))
    return bwt, indices

# wrapper for the processing functions used to compute
# auxiliary data structures for efficient BWT matching
# see file `preprocess_bwt.py`
def _preprocess_bwt(bwt):
    first_occurence = _get_first_occurence_fn(bwt)
    count = _get_count_fn(bwt)
    return first_occurence, count

# class encapsulating exact matching with Burrows-Wheeler transform
#
# Fields:
#   _text: string, the target string
#   _bwt: string, the burrows-wheeler transform of target string
#   _suffix_array: [int], suffix array of target string
#   first_occurence: function returning first occurence of each symbol in
#                     first column of sorted rotation table for bwt, see below
#   count: function returning number of occurences of each symbol up to
#           a given position, see below
#
# Notes:
#   After initializing: `bwt = BWT(target)`:
#
#   `bwt.first_occurence(symbol)` returns the row in which symbol occurs first
#       in the first column of the sorted rotation table corresponding to the BWT
#       of target string
#
#   `bwt.count(symbol, position)` returns the number of occurrences of symbol
#       up to given position in BWT of target string
class BWT:
    def __init__(self, target=None):
        self._text = target
        self._bwt = None
        self._suffix_array = None
        if target is not None:
            self._bwt, self._suffix_array = _construct(self._text)
        self.first_occurence = None
        self.count = None
        if target is not None:
            self.first_occurence, self.count = _preprocess_bwt(self._bwt)

    # THIS IS A STUB, YOU NEED TO IMPLEMENT THIS
    #
    # return indices for positions in target string that match
    # query exactly
    #
    # Input:
    #   pattern: string, query string
    #
    # Output:
    #   [int], array of indices of exact matches of query in target
    #          array is empty if no exact matches found
    def get_matches(self, pattern):
        top, bottom = self._get_matching_rows(pattern)

        if top == -1:
            return []

        # return matching indices based on top and bottom pointers
        # YOU NEED TO FILL THIS IN
        return self._suffix_array[top:(bottom+1)]

    # find pointer in first column of last occurence of symbol
    # before 'index' in last column
    def _move_pointer(self, symbol, index):
        return self.first_occurence(symbol) + self.count(symbol, index)

    # THIS IS A STUB, YOU NEED TO IMPLEMENT THIS
    #
    # return top, bottom pointers for rows of sorted rotations table
    # that start with query
    #
    # Input:
    #   pattern: string, query string
    #
    # Output:
    #   tuple (top, bottom): top and bottom pointers for consecutive rows in
    #       sorted rotations table that start with exact matches to query string
    #       returns (-1, -1) if no matches are found
    def _get_matching_rows(self, pattern):
        # have top and bottom pointers at first and last
        # row of M matrix
        top = 0
        bottom = len(self._bwt) - 1

        # while there are rows to be checked
        while top <= bottom:
            # and patterns symbols to be checked
            if len(pattern) > 0:
                symbol = pattern[-1] # check the last symbol in pattern
                pattern = pattern[:-1] # remove last symbol in pattern

                # update top and bottom pointers
                top = self._move_pointer(symbol, top)
                bottom = self._move_pointer(symbol, bottom + 1) - 1
            else:
                return (top,bottom)
        return (-1, -1)
