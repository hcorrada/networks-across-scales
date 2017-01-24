import sys

# compute the probability of kmer given profile matrix
def profile_probability(kmer, profile):
    prob = 1.
    k = len(kmer)
    for i in xrange(k):
        # probability of k-mer's nucleotide at position i
        # is given by 'kmer[i]' row of the ith column of matrix
        prob *= profile[kmer[i]][i]
    return prob

# find the most probable kmer in text given profile matrix
def most_probable_kmer(text, k, profile):
    # assume first kmer is most probable
    best_kmer = text[:k]
    best_probability = profile_probability(best_kmer, profile)

    # check the remaining kmers
    n = len(text)
    for i in xrange(1, n-k+1):
        kmer = text[slice(i,i+k)]
        prob = profile_probability(kmer, profile)
        if prob >= best_probability:
            best_probability = prob
            best_kmer = kmer
    return best_kmer

filename = sys.argv[1]
with open(filename, 'r') as f:
    text = f.readline().strip()
    k = int(f.readline().strip())
    profile = {}
    indices = list('ACGT')
    for i in xrange(4):
        string = f.readline().strip()
        profile[indices[i]] = map(float, string.split())

    print most_probable_kmer(text, k, profile)
