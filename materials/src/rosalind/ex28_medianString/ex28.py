import sys
import itertools

# compute edit distance between two strings
# input:
#  pattern: a string
#  string: another string
# output:
#  the edit distance between the two strings
def edit_distance(pattern, string):
    count = 0
    for i in xrange(len(pattern)):
        count += pattern[i] != string[i]
    return count

# find the k-mer in string with smallest edit distance to pattern
# input:
#  pattern: a k-mer
#  string: a longer string in which to search for minimum distance k-mer
# output:
#  the smallest edit distance 
def min_distance(pattern, string):
    k = len(pattern)
    best_distance = k+1
    n = len(string)
    
    for i in xrange(n-k+1):
        d = edit_distance(pattern, string[slice(i, i+k)])
        if d <= best_distance:
            best_distance = d
    return best_distance

# compute the distance between a k-mer and *a set* of strings
# defined as the sum of the minimum edit distance between the given
# k-mer, and a k-mer in each of the strings
# input:
#  pattern: a k-mer
#  strings: a set of strings
# output:
#   the sum of minimum distances
def dist(pattern, strings):
    dist = 0
    for string in strings:
        dist += min_distance(pattern, string)
    return dist

# return the median k-mer for a given set of strings
# median k-mer is defined as the k-mer that minimizes it's distance to the set of strings
def median_string(k, strings):
    t = len(strings)
    best_pattern = ''
    best_distance = k * t

    # generate all possible k-mers
    for pattern in itertools.product('ACGT', repeat=k):
        pattern = ''.join(pattern)
        # compute the distance of pattern
        d = dist(pattern, strings)

        # check if it's the best so far
        if d <= best_distance:
            best_distance = d
            best_pattern = pattern
    return best_pattern

filename = sys.argv[1]
with open(filename, 'r') as f:
    k = int(f.readline().strip())
    strings = []
    for string in f:
        strings.append(string.strip())
    res = median_string(k, strings)
    print res
