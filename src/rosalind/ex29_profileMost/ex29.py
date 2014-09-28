import sys

def profile_probability(kmer, profile):
    prob = 1.
    k = len(kmer)
    for i in xrange(k):
        prob *= profile[kmer[i]][i]
    return prob

def most_probable_kmer(text, k, profile):
    best_kmer = ''
    best_probability = 0.

    n = len(text)
    for i in xrange(n-k+1):
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
    
