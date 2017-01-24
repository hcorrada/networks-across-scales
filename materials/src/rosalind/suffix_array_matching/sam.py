import sys
import math

# construct suffix array from given text
def construct_suffix_array(text):
    # setup list of indices
    sa = range(len(text))

    # sort the indices by corresponding suffixes
    sa = sorted(sa, key=lambda i: text[i:])
    return sa

def find_matches(pattern, sa, text):
    min_index = 0
    max_index = len(text) - 1

    while min_index < max_index:
        mid_index = (min_index + max_index) / 2
        cmp_result = ">" if pattern > text[sa[mid_index]:] else "<="
        #print "%d,%d,%d" % (min_index, mid_index, max_index),
        #print "%s:%s:%s" % (pattern, text[sa[mid_index]:], cmp_result)

        if pattern > text[sa[mid_index]:]:
            min_index = mid_index + 1
        else:
            max_index = mid_index
    first = min_index
    #print "first: %d:%d, %s:%s" % (first, sa[first], pattern, text[sa[first]:])

    min_index = first + 1
    max_index = len(text) - 1
    while min_index < max_index:
        mid_index = (min_index + max_index) / 2
        cmp_result = ">" if pattern > text[sa[mid_index]:] else "<="
        #print "%d,%d,%d" % (min_index, mid_index, max_index),
        #print "%s:%s:%s" % (pattern, text[sa[mid_index]:], cmp_result)

        if pattern < text[sa[mid_index]:]:
            max_index = mid_index
        else:
            min_index = mid_index + 1
    last = max_index
    #print "last: %d:%d, %s:%s" % (last, sa[last], pattern, text[sa[last]:])

    res = range(first, last + 1) if first < last else []
    return res

def readdat(filename):
    with open(filename, 'r') as f:
        all_lines = f.readlines()
        text = all_lines[0].strip()
        patterns = map(lambda s: s.strip(), all_lines[1:])
    return text, patterns

def naive_find(text, pattern):
    res = []
    for i in xrange(len(text) - len(pattern) + 1):
        if pattern == text[slice(i,i+len(pattern))]:
            res.append(i)
    return res

def main(filename):
    text, patterns = readdat(filename)
    text += "$"
    sa = construct_suffix_array(text)

    res = []
    for pattern in patterns:
        matches = find_matches(pattern, sa, text)
        matches = map(lambda i: sa[i], matches)
        matches2 = sorted(naive_find(text, pattern))
        print matches, matches2
        print map(lambda i: text[i:], matches)
        print map(lambda i: text[slice(i,i+len(pattern))], matches2)

        res += matches
    print " ".join(map(str, sorted(res)))

if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
