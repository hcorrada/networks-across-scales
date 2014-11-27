# construct suffix array from given text
def make_suffix_array(text):
    # setup list of indices
    sa = [i for i in xrange(len(text))]

    # sort the indices by corresponding suffixes
    sa = sorted(sa, key=lambda i: text[i:])
    return sa

def make_partial_suffix_array(text, k):
    sa = make_suffix_array(text)
    sa = zip([i for i in xrange(len(sa))], sa)
    sa = filter(lambda x: x[1] % k == 0, sa)
    return sa

    
