# construct suffix array from given text
def make_suffix_array(text):
    # setup list of indices
    sa = [i for i in xrange(len(text))]

    # sort the indices by corresponding suffixes
    sa = sorted(sa, key=lambda i: text[i:])
    return sa

# class encapsulating suffix array
class SuffixArray:
    def __init__(self, text):
        self._sa = make_suffix_array(text)

    # given index into suffix array, return
    # index of suffix in original text
    def get_match_index(self, index):
        return self._sa[index]

    # return suffix indices for given range
    # of suffix array indices
    def get_match_indices(self, top, bottom):
        out = []
        for i in xrange(top,bottom+1):
            out.append(self.get_match_index(i))
        return out

# construct a partial suffix array
# only storing entries in suffix array that are
# multiples of k
def make_partial_suffix_array(text, k):
    # make the standard suffix array
    sa = make_suffix_array(text)

    # keep only entries that are multiples of k
    sa = zip([i for i in xrange(len(sa))], sa)
    sa = filter(lambda x: x[1] % k == 0, sa)

    # return indices kept and values as lists
    indices, values = zip(*sa)
    return list(indices), list(values)

# class encapsulating a partial suffix array
# needs a bwt object to compute index of arbitrary suffixes
class PartialSuffixArray(SuffixArray):
    def __init__(self, text, k, bwt):
        self._indices, self._values = make_partial_suffix_array(text, k)
        self._bwt = bwt

    # given index into full suffix array,
    # return index of suffix in original text
    # needs bwt to do this using partial suffix array
    def get_match_index(self, index):
        # do rotations until an index in partial suffix
        # array is found
        steps = 0
        while index not in self._indices:
            steps += 1
            index = self._bwt.move_back(index)
        index = self._indices.index(index)
        return steps + self._values[index]


    
