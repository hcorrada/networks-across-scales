import sys

# create a function that computes the edit distance
# to given pattern
# input:
#  pattern: a string
#
# output:
#  a function that computes edit distance to patter
def edit_distance(pattern):
    pattern_list = list(pattern)
    return lambda s: sum([x != y for (x,y) in zip(pattern_list, list(s))])

# approximate pattern matching
# input:
#  pattern: the pattern to match
#  text: the strings to match pattern to
#  d: maximum number of mismatches 
def approximate_pattern_match(pattern, text, d):
    k = len(pattern)
    n = len(text)

    # create a function that computes edit distance
    # to given pattern
    edfun = edit_distance(pattern)

    indices = []

    # check if position in text
    for i in xrange(n-k+1):
        # add to list of matches if edit distances is less than or equal
        # to maximum number of allowed mismatches
        if edfun(text[slice(i,i+k)]) <= d:
            indices.append(i)
            
    return indices

filename = sys.argv[1]
with open(filename, 'r') as f:
    pattern = f.readline().strip()
    text = f.readline().strip()
    d = int(f.readline().strip())
    print " ".join(map(str, approximate_pattern_match(pattern, text, d)))
    
            
