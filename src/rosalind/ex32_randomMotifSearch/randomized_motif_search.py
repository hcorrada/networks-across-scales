import sys
import random
import math

# utility structure
alphabet = list('ACGT')

# compute the entropy of a profile
def profile_entropy(profile, k):
    score = 0
    for j in xrange(k):
        cur_ent = 0
        for nuc in alphabet:
            cur_ent += -( profile[nuc][j] * math.log(profile[nuc][j], 2) )
        score += cur_ent
    return score

# construct a profile using pseudocounts
# 'indices': a list of starting indices for k-mers used to estimate profile
def make_profile(indices, dna, k):
    t = len(indices)

    # initialize profile
    profile = {}
    for nuc in alphabet:
        profile[nuc] = []

    # for each column in profile
    for i in xrange(k):
        # initialize counts with pseudocount = 1
        counts = dict.fromkeys(alphabet, 1.)

        # add count for each kmer
        for j in xrange(t):
            nuc = dna[j][indices[j]+i]
            counts[nuc] += 1.

        # compute probability for each nucleotide
        for nuc, count in counts.iteritems():
            profile[nuc].append( count / ( t + len(alphabet)) )
    return profile

# compute the probability of given k-mer
# based on given profile
def profile_probability(kmer, profile):
    prob = 1.
    k = len(kmer)
    for i in xrange(k):
        prob *= profile[kmer[i]][i]
    return prob

# find the index of most probable k-mer in text
# based on given profile
def most_probable_index(text, profile, k):
    best_index = None
    best_probability = -1.

    n = len(text)
    for i in xrange(n-k+1):
        # get kmer at current position
        kmer = text[slice(i, i+k)]

        # compute probability of kmer
        prob = profile_probability(kmer, profile)

        # update index if needed
        if prob > best_probability:
            best_probability = prob
            best_index = i
    return best_index

# get indices of most-probable kmers in t
# based on profile
def get_new_indices(profile, dna, t, k):
    # initialize all indices to 0
    indices = [0] * t

    # get index of most probable kmer for each
    # string in 'dna'
    for i in xrange(t):
        indices[i] = most_probable_index(dna[i], profile, k)
    return indices

# get kmers starting at each index from set of
# strings 'dna'
def get_kmers(indices, dna, k, t):
    kmers = [""] * t

    for i in xrange(t):
        start_index = indices[i]
        kmers[i] = dna[i][slice(start_index, start_index + k)]
    return kmers

# find one candidate profile
# using randomized search
def find_one_profile(k, t, dna):
    n = len(dna[0])

    # initialize to set of random indices
    indices = [random.randrange(n-k+1) for _ in xrange(t)]
    best_score = 1. * k * t

    # check if current set of indices can
    # be improved
    while True:
        # construct a profile given current set of indices
        profile = make_profile(indices, dna, k)

        # compute the entropy of profile
        score = profile_entropy(profile, k)

        # update minimum entropy if needed
        if score < best_score:
            best_score = score
            indices = get_new_indices(profile, dna, t, k)
        else:
            # need this here since it was overwritten
            # before scoring
            profile = make_profile(indices, dna, k)
            return best_score, profile, indices

# outer loop for random motif search
# performs random search for 1000 random choices
def random_motif_search(k, t, dna):
    n = 1000
    best_score = 1. * k * t
    best_indices = None

    for _ in xrange(n):
        cur_score, _, cur_indices = find_one_profile(k, t, dna)

        if cur_score < best_score:
            best_indices = cur_indices
            best_score = cur_score
    return get_kmers(best_indices, dna, k, t)

# read input and call search function
filename = sys.argv[1]
with open(filename, 'r') as f:
    k, t = map(int, f.readline().strip().split())
    dna = []
    for i in xrange(t):
        dna.append(f.readline().strip())
    print '\n'.join(random_motif_search(k, t, dna))
