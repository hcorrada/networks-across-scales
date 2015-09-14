import sys

# generate the d neighborhood of a DNA string
nucleotides = list('ACGT')

def neighbors(pattern, d):
    if d == 0:
        return [(pattern,0)]
    if len(pattern) == 1:
        return [(nuc, 1 if nuc != pattern else 0) for nuc in nucleotides]

    neighborhood = set()
    suffix_neighbors = neighbors(pattern[1:], d)
    for (text, cur_d) in suffix_neighbors:
        if cur_d < d:
            for nuc in nucleotides:
                distance = cur_d + 1 if nuc != pattern[0] else cur_d
                neighborhood.add((nuc + text, distance))
        else:
            neighborhood.add((pattern[0] + text, d))
    return neighborhood

def read_input(f):
    pattern = f.readline().strip()
    d = int(f.readline().strip())
    return pattern, d

filename = sys.argv[1]
with open(filename, 'r') as f:
    pattern, d = read_input(f)
    neighborhood = neighbors(pattern, d)
    for neighbor, dist in neighborhood:
        print neighbor
