import sys
from collections import defaultdict

def count_kmers(text, k):
    counts = defaultdict(int)
    n = len(text)
    for i in xrange(n-k+1):
        kmer = text[slice(i, i+k)]
        counts[kmer] += 1
    return counts

def clump_finder(genome, k, L, t):
    clumps = set()
    n = len(genome)

    counts = count_kmers(genome[slice(0,L)], k)
    for kmer, kmer_count in counts.items():
        if kmer_count >= t:
            clumps.add(kmer)

    for i in xrange(1, n - L + 1):
        first_kmer = genome[slice(i-1, i-1+k)]
        counts[first_kmer] -= 1
        last_kmer = genome[slice(i+L-k, i+L)]
        counts[last_kmer] += 1
        if counts[last_kmer] >= t:
            clumps.add(last_kmer)
    return clumps

filename = sys.argv[1]
with open(filename, 'r') as f:
    genome = f.readline().strip()
    k, L, t = map(int, f.readline().strip().split())
    res = clump_finder(genome, k, L, t)
    print ' '.join(res)
