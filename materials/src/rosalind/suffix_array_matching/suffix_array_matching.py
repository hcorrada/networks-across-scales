import sys
# make a comparison function based on suffixes of text
def make_cmp_fn(text):
    # comparison function, takes two suffix indices
    def mycmp(x,y):
        # how many positions to check (length of text - largest index)
        n = len(text) - max(x,y)

        # start at the first position of each suffix
        d = cmp(text[x], text[y])
        i = 1

        # keep going while suffixes are equal
        while d == 0 and i < n:
            # compare characters at i'th position
            d = cmp(text[x+i], text[y+i])
            i += 1
        # return the last comparison
        return d
    # return the comparison function
    return mycmp

# construct suffix array from given text
def construct_suffix_array(text):

    # setup list of indices
    sa = range(len(text))

    # sort the indices by corresponding suffixes
    sa = sorted(sa, cmp=make_cmp_fn(text))
    return sa

def make_pattern_cmp_fn(text):
    def mycmp(index, pattern):
        n = min(len(text) - index, len(pattern))
        d = cmp(text[index], pattern[0])
        i = 1
        while d == 0 and i < n:
            d = cmp(text[index+i], pattern[i])
            i += 1
        if d == 0 and i == len(pattern):
            return 1
        return d
    return mycmp

def find_matches(pattern, sa, text):
    cmp = make_pattern_cmp_fn(text)
    min_index = 0
    max_index = len(text) - 1

    while min_index < max_index:
        mid_index = (min_index + max_index) / 2
        cmp_result = cmp(sa[mid_index], pattern)
        print "%d,%d,%d" % (min_index, mid_index, max_index),
        print "%s:%s:%d" % (pattern, text[sa[mid_index]:], cmp_result)
        print
        if cmp_result < 0:
            min_index = mid_index + 1
        else:
            max_index = mid_index
    first = min_index
    print "first: %d, %s:%s" % (first, pattern, text[sa[first]:])
    max_index = len(text) - 1

    while min_index < max_index:
        mid_index = (min_index + max_index) / 2
        cmp_result = cmp(sa[mid_index], pattern)
        print "%d,%d,%d" % (min_index, mid_index, max_index),
        print "%s:%s:%d" % (pattern, text[sa[mid_index]:], cmp_result)
        print

        if cmp_result > 0:
            max_index = mid_index
        else:
            min_index = mid_index + 1
    last = max_index
    print "last: %d, %s:%s" % (first, pattern, text[sa[first]:])

    result = [] if first > last else sa[first:last+1]
    return result

def readdat(filename):
    with open(filename, 'r') as f:
        all_lines = f.readlines()
        text = all_lines[0].strip()
        patterns = all_lines[1:]
    return text, patterns

def main(filename):
    text, patterns = readdat(filename)
    text += "$"
    sa = construct_suffix_array(text)

    res = []
    for pattern in patterns:
        pattern = pattern.strip()
        matches = find_matches(pattern, sa, text)
        res += matches
    print " ".join(map(str, sorted(res)))

if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
