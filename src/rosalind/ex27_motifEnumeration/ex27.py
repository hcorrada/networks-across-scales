import sys
import itertools

def find_kmers(text, k):
    res = []
    n = len(text)
    for i in xrange(n-k+1):
        kmer = text[slice(i, i+k)]
        if kmer not in res:
            res.append(kmer)
    return res

def build_neighborhood(pattern, d):
    res = []

    for i in itertools.combinations(range(len(pattern)), d):
        for j in itertools.product('ACGT', repeat=d):
            temp=list(pattern)
            for k in xrange(d):
                temp[i[k]] = j[k]
            temp = "".join(temp)
            if temp not in res:
                res.append(temp)
    return res

def motif_enumeration(k, d, strings):
    final_tab = {}
    for string in strings:
        cur_tab = []
        kmers = find_kmers(string, k)
        for kmer in kmers:
            neighbors = build_neighborhood(kmer, d)
            for candidate in neighbors:
                if candidate not in cur_tab:
                    cur_tab.append(candidate)
                    if candidate in final_tab:
                        final_tab[candidate] += 1
                    else:
                        final_tab[candidate] = 1

    res = []
    for kmer, count in final_tab.iteritems():
        if count == len(strings):
            res.append(kmer)
    return res

filename = sys.argv[1]
with open(filename, 'r') as f:
    k, d = map(int, f.readline().strip().split())
    strings = []
    for string in f:
        strings.append(string.strip())
    print " ".join(motif_enumeration(k, d, strings))
    
