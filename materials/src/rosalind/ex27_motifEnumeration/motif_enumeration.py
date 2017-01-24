import sys

# generate the d neighborhood of a DNA string
nucleotides = list('ACGT')

def get_neighbors(pattern, d):
    if d == 0:
        return [(pattern,0)]
    if len(pattern) == 1:
        return [(nuc, 1 if nuc != pattern else 0) for nuc in nucleotides]

    neighborhood = set()
    suffix_neighbors = get_neighbors(pattern[1:], d)
    for (text, cur_d) in suffix_neighbors:
        if cur_d < d:
            for nuc in nucleotides:
                distance = cur_d + 1 if nuc != pattern[0] else cur_d
                neighborhood.add((nuc + text, distance))
        else:
            neighborhood.add((pattern[0] + text, d))
    return neighborhood

def hamming_distance(pattern, target):
    d = 0
    k = len(pattern)
    for i in xrange(k):
        if pattern[i] != target[i]:
            d += 1
    return d

# check if pattern occurs at least once in string
# with at most d mismatches
def is_implanted(pattern, string, d):
    k = len(pattern)
    n = len(string)
    for i in xrange(n-k+1):
        # check if distance between pattern and k-mer at
        # position i is less than d
        if hamming_distance(pattern, string[slice(i,i+k)]) <= d:
            # it is, found an occurence
            return True
    # didn't find any occurence within distance d in this string
    return False

# enumerate all (k,d) motifs in DNA strings
# input:
#   k: motif length
#   d: maximum mismatches
#   strings: list of DNA strings
# output:
#   list of k-mers that occur at least once in
#     all strings with at most d mismatches
def motif_enumeration(k, d, strings):
    patterns = set()

    for string in strings:
        n = len(string)
        for i in xrange(n-k+1):
            kmer = string[slice(i,i+k)]
            neighbors = get_neighbors(kmer, d)
            for neighbor, _ in neighbors:
                if neighbor in patterns:
                    continue
                # check if neighbor is implanted on all strings
                if all(map(lambda x: is_implanted(neighbor, x, d), strings)):
                    patterns.add(neighbor)
    return patterns

filename = sys.argv[1]
with open(filename, 'r') as f:
    k, d = map(int, f.readline().strip().split())
    strings = []
    for string in f:
        strings.append(string.strip())
    print " ".join(motif_enumeration(k, d, strings))
