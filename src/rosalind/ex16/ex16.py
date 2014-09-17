import sys

def edit_distance(pattern):
    pattern_list = list(pattern)
    return lambda s: sum([x != y for (x,y) in zip(pattern_list, list(s))])

def approximate_pattern_match(pattern, text, d):
    k = len(pattern)
    n = len(text)
    edfun = edit_distance(pattern)

    indices = []
    
    for i in xrange(n-k+1):
        if edfun(text[slice(i,i+k)]) <= d:
            indices.append(i)
            
    return indices

filename = sys.argv[1]
with open(filename, 'r') as f:
    pattern = f.readline().strip()
    text = f.readline().strip()
    d = int(f.readline().strip())
    print " ".join(map(str, approximate_pattern_match(pattern, text, d)))
    
            
