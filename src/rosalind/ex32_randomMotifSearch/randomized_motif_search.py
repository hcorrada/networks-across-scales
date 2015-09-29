import sys
import random
import math

alphabet = list('ACGT')

def profile_score (profile, k):
    score = 0
    for j in xrange(k):
        cur_max = -1.
        for nuc in alphabet:
            p = profile[nuc][j]
            if p > cur_max:
                cur_max = p
        score += 1-p

def profile_entropy(profile, k):
    score = 0
    for j in xrange(k):
        cur_ent = 0
        for nuc in alphabet:
            cur_ent += -( profile[nuc][j] * math.log(profile[nuc][j], 2) )
        score += cur_ent
    return score

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

def profile_probability(kmer, profile):
    prob = 1.
    k = len(kmer)
    for i in xrange(k):
        prob *= profile[kmer[i]][i]
    return prob

def most_probable_index(text, profile, k):
    best_index = None
    best_probability = -1.

    n = len(text)
    for i in xrange(n-k+1):
        kmer = text[slice(i, i+k)]
        prob = profile_probability(kmer, profile)
        if prob > best_probability:
            best_probability = prob
            best_index = i
    return best_index

def get_new_indices(profile, dna, t, k):
    indices = [0] * t
    for i in xrange(t):
        indices[i] = most_probable_index(dna[i], profile, k)
    return indices

def get_kmers(indices, dna, k, t):
    kmers = [""] * t

    for i in xrange(t):
        start_index = indices[i]
        kmers[i] = dna[i][slice(start_index, start_index + k)]
    return kmers

def find_one_profile(k, t, dna):
    n = len(dna[0])
    indices = [random.randrange(n-k+1) for _ in xrange(t)]
    best_score = 1. * k * t

    while True:
        profile = make_profile(indices, dna, k)
        score = profile_entropy(profile, k)
        if score < best_score:
            best_score = score
            indices = get_new_indices(profile, dna, t, k)
        else:
            return best_score, profile, indices

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

filename = sys.argv[1]
with open(filename, 'r') as f:
    k, t = map(int, f.readline().strip().split())
    dna = []
    for i in xrange(t):
        dna.append(f.readline().strip())
    print '\n'.join(random_motif_search(k, t, dna))
