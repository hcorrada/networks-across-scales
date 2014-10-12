import sys
import itertools

def edit_distance(pattern, string):
    count = 0
    for i in xrange(len(pattern)):
        count += pattern[i] != string[i]
    return count

def min_distance(pattern, string):
    k = len(pattern)
    best_distance = k+1
    n = len(string)
    
    for i in xrange(n-k+1):
        d = edit_distance(pattern, string[slice(i, i+k)])
        if d <= best_distance:
            best_distance = d
    return best_distance

def dist(pattern, strings):
    dist = 0
    for string in strings:
        dist += min_distance(pattern, string)
    return dist

def median_string(k, strings):
    t = len(strings)
    best_pattern = ''
    best_distance = k * t

    for pattern in itertools.product('ACGT', repeat=k):
        pattern = ''.join(pattern)
        d = dist(pattern, strings)

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
