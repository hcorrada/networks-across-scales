# construct suffix array from given text
def make_suffix_array(text):
    # setup list of indices
    sa = [i for i in xrange(len(text))]

    # sort the indices by corresponding suffixes
    sa = sorted(sa, key=lambda i: text[i:])
    return sa

class SuffixArray:
    def __init__(self, text):
        self._sa = make_suffix_array(text)

    def get_match_index(self, index):
        return self._sa[index]
    
    def get_match_indices(self, top, bottom):
        out = []
        for i in xrange(top,bottom+1):
            out.append(self.get_match_index(i))
        return out
    
def make_partial_suffix_array(text, k):
    sa = make_suffix_array(text)
    sa = zip([i for i in xrange(len(sa))], sa)
    sa = filter(lambda x: x[1] % k == 0, sa)
    indices, values = zip(*sa)
    return list(indices), list(values)

class PartialSuffixArray(SuffixArray):
    def __init__(self, text, k, bwt):
        self._indices, self._values = make_partial_suffix_array(text, k)
        self._bwt = bwt

    def get_match_index(self, index):
        steps = 0
        while index not in sa_indices:
            steps += 1
            symbol = bwt[index]
            index = first_occurence[symbol] + count(symbol, index+1)
        return steps + sa_values[index]


    
