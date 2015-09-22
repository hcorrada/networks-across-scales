import sys
from collections import defaultdict

# count k-mers in given text
# input:
#   text: a string (of DNA)
#   k: k-mer length
# output:
#   a dictionary index by k-mer, containg each k-mers count
def count_kmers(text, k):
    counts = defaultdict(int)
    n = len(text)
    for i in xrange(n-k+1):
        kmer = text[slice(i, i+k)]
        counts[kmer] += 1
    return counts

# find kmers with occurences that clump: occurs at least t times in a window of length L
# input:
#   genome: string of DNA
#   k: k-mer length
#   L: window length
#   t: number of occurences
def clump_finder(genome, k, L, t):
    clumps = set()
    n = len(genome)

    # count kmers in the first L-long window
    # add add to clump set if it occurs more than t times
    counts = count_kmers(genome[slice(0,L)], k)
    for kmer, kmer_count in counts.items():
        if kmer_count >= t:
            clumps.add(kmer)

    # now check all remaining windows
    for i in xrange(1, n - L + 1):
        # find the first k-mer in the window
        # and adjust its count (1 fewer occurence)
        first_kmer = genome[slice(i-1, i-1+k)]
        counts[first_kmer] -= 1

        # find the last k-mer in the window
        # and adjust its count (1 more occurrence)
        last_kmer = genome[slice(i+L-k, i+L)]
        counts[last_kmer] += 1

        # check if the last_kmer clumps
        if counts[last_kmer] >= t:
            clumps.add(last_kmer)
    return clumps

filename = sys.argv[1]
with open(filename, 'r') as f:
    genome = f.readline().strip()
    k, L, t = map(int, f.readline().strip().split())
    res = clump_finder(genome, k, L, t)
    print ' '.join(res)
