import sys
import random
import math

# seed random number generator
random.seed(2)

# list of nucleotides
alphabet = list('ACGT')

def profile_score(profile, k):
    score = 0.
    for j in xrange(k):
        curmax = -1.
        for nuc in alphabet:
            vec = profile[nuc]
            if vec[j] > curmax:
                curmax = vec[j]
        score += 1. - curmax

# compute the entropy of a profile
def profile_entropy(profile, k):
    score = 0.
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

# select index into list with probability proportional to weight
# input:
#  weights: a list of weights
# output:
#  a random selection from list
def random_choice(weights):
    # normalize weights so they sum to 1
    total_weight = 1. * sum(weights)
    nitems = len(weights)
    probs = []

    for i in xrange(nitems):
        probs.append(weights[i] / total_weight)

    # generate random number between 0 and 1
    u = random.uniform(0, 1.)

    # find the chosen index
    for i in xrange(nitems):
        if u <= probs[i]:
            break
        u -= probs[i]
    return i


# compute the probability of given k-mer
# based on given profile
def profile_probability(kmer, profile):
    prob = 1.
    k = len(kmer)
    for i in xrange(k):
        prob *= profile[kmer[i]][i]
    return prob

# sample a starting position from text, weighted by profile probability
# input:
#  text: a string
#  profile: a profile
#  k: k-mer length
# output
#  a randomly chosen position within text, with probability given by profile
def profile_sample(text, profile, k):
    n = len(text)

    # compute profile probability for each k-mer
    probs = [profile_probability(text[slice(i, i+k)], profile) for i in xrange(n-k+1)]

    # randomly choose a starting position with probability proportional to probs
    index = random_choice(probs)

    # return chosen k-mer index
    return index

# one run of the gibbs sampler
# input:
#  k: k-mer length
#  t: number of strings in dna
#  N: number of gibbs sampler moves to try
#  dna: list of strings
# output:
#  score: entropy of the best profile found
#  profile: the best profile found
#  indices: the starting positions of k-mers used to estimate profile
def find_one_profile(k, t, N, dna):
    n = len(dna[0])

    # generate random index vector
    indices = [random.randrange(n-k+1) for _ in xrange(t)]

    # and score it
    profile = make_profile(indices, dna, k)
    best_score = profile_entropy(profile, k)
    best_indices = indices

    old_index = None
    string_to_change = None

    # try N steps
    for i in xrange(N):
        # choose a string in dna to resample
        string_to_change = random.randrange(t)

        # make profile from current k-mers except chosen string
        cur_indices = indices[:string_to_change] + indices[(string_to_change+1):]
        cur_dna = dna[:string_to_change] + dna[(string_to_change+1):]
        cur_profile = make_profile(cur_indices, cur_dna, k)

        # now sample a new index with probability given by current profile
        new_index = profile_sample(dna[string_to_change], cur_profile, k)
        old_index = indices[string_to_change]
        indices = indices[:string_to_change] + [new_index] + indices[(string_to_change+1):]

        new_profile = make_profile(indices, dna, k)
        new_score = profile_entropy(new_profile, k)

        if new_score < best_score:
            best_indices = indices
            best_score = new_score
            profile = new_profile
        else:
            indices[string_to_change] = old_index

    return best_score, profile, best_indices

# get kmers starting at each index from set of
# strings 'dna'
def get_kmers(indices, dna, k, t):
    kmers = [""] * t

    for i in xrange(t):
        start_index = indices[i]
        kmers[i] = dna[i][slice(start_index, start_index + k)]
    return kmers

# gibbs sampler with multiple restarts
#  k: k-mer length
#  t: number of strings in dna
#  N: number of gibbs sampler moves to try
#  dna: list of strings
# output:
#  a set lof k-mers from dna
def gibbs_sampling_search(k, t, N, dna):
    # number of restarts to try
    num_starts = 20

    # initialize score
    best_score = 1. * k * t
    best_indices = None

    # run gibbs sampler multiple times
    for _ in xrange(num_starts):
        cur_score, _, cur_indices = find_one_profile(k, t, N, dna)

        # keep score for new motifs if they improve the score
        if cur_score < best_score:
            best_indices = cur_indices
            best_score = cur_score
    return get_kmers(best_indices, dna, k, t)

filename = sys.argv[1]
with open(filename, 'r') as f:
    k, t, N = map(int, f.readline().strip().split())
    dna = []
    for i in xrange(t):
        dna.append(f.readline().strip())
    print '\n'.join(gibbs_sampling_search(k, t, N, dna))
