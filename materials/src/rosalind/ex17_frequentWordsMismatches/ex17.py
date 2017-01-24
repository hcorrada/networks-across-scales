import sys
import itertools

def count_kmers(text, k):
    n = len(text)
    tab = {}
    for i in xrange(n-k+1):
        kmer = text[slice(i, i+k)]
        tab[kmer] = tab[kmer] + 1 if kmer in tab else 1
    return tab

def build_neighborhood(pattern, d):
    res = {}
    k = range(len(pattern))
    for i in itertools.combinations(k, d):
        for j in itertools.product('ACGT', repeat=d):
            temp=list(pattern)
            for k in xrange(d):
                temp[i[k]] = j[k]
            temp = "".join(temp)
            if temp not in res:
                res[temp] = 0
    return res.keys()

def process_kmer(pattern, final_tab, kmer_counts, d):
    neighbors = build_neighborhood(pattern, d)
    for kmer in neighbors:
        if kmer in final_tab:
            final_tab[kmer] += kmer_counts[pattern]
        else:
            final_tab[kmer] = kmer_counts[pattern]
    return final_tab

def freq_words_mismatch(text, k, d):
    kmer_counts = count_kmers(text, k)
    kmers = kmer_counts.keys()

    final_tab = {}
    for kmer in kmers:
        final_tab = process_kmer(kmer, final_tab, kmer_counts, d)

    maxcount = 0
    res = []

    for (kmer, count) in final_tab.iteritems():
        if count > maxcount:
            maxcount = count
            res = [kmer]
        elif count == maxcount:
            res.append(kmer)
    return res
    
filename = sys.argv[1]
with open(filename, 'r') as f:
    text = f.readline().strip()
    k, d = map(int, f.readline().strip().split())
    print " ".join(freq_words_mismatch(text, k, d))
    
