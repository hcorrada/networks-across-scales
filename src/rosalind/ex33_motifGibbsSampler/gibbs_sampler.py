import sys
import random

# seed random number generator 
random.seed(2)

# list of nucleotides
alphabet = list('ACGT')

# score a set of k-mers
# input:
#  a set of k-mers
# output:
#  the 'conservation' score
def motif_score(motifs):
    k = len(motifs[0])
    t = len(motifs)
    
    score = 0
    # for each position in motif (column in table of strings)
    for i in xrange(k):
        # count the number of times each nucleotide occurs in this position
        # across given strings
        count = dict.fromkeys(alphabet, 0)
        for j in xrange(t):
            count[motifs[j][i]] += 1
        # find the most common (consensus) nucleotide        
        most_common = max(count.items(), key = lambda x: x[1])[0]
        # add the number of nucleotides that are not equal to the consensus nucleotide to the score        
        cur_score = sum(count.values()) - count[most_common]
        score += cur_score
    return score

# construct a profile given a set of k-mers
# input:
#  motifs: a set of k-mers
#  k: k-mer length
# output:
#  a profile: hash table, keys are nucleotides, values are list of frequencies
def make_profile(motifs, k):
    t = len(motifs)

    # initialize profile
    profile = {}
    for nuc in alphabet:
        profile[nuc] = []

    # for each column in profile
    for i in xrange(k):
        # initialize counts with pseudocount = 1
        counts = dict.fromkeys(alphabet, 1.)

        # add count for each motif
        for j in xrange(t):
            counts[motifs[j][i]] += 1.

        # compute probability for each nucleotide
        for nuc, count in counts.iteritems():
            profile[nuc].append( count / ( t + len(alphabet)) )
    return profile

# select item from list with probability proportional to weight
# input:
#  weights: a list of weights
# output:
#  a random selection from list
def random_choice(weights):
    # normalize weights so they sum to 1
    total_weight = sum(weights)
    nitems = len(weights)
    for i in xrange(nitems):
        weights[i] /= total_weight

    # generate random number between 0 and 1
    choice = -1
    u = random.uniform(0, 1.)

    # find the chosen item
    while u >= 0:
        choice += 1
        u -= weights[choice]
    return choice

# compute the probabilty of a given k-mer according to given profile
# input:
#  kmer: a k-mer
#  profile: a profile
# output:
#  the probability
def profile_probability(kmer, profile):
    prob = 1.
    k = len(kmer)
    for i in xrange(k):
        # for each character in kmer, multiply in the probability of that nucleotide
        # according to profile        
        prob *= profile[kmer[i]][i]
    return prob

# sample a starting position from text, weighted by profile probability
# input:
#  text: a string
#  profile: a profile
#  k: k-mer length
# output
#  a randomly chosen k-mer from text, with probability given by profile
def profile_sample(text, profile, k):
    n = len(text)
    # compute profile probability for each k-mer
    probs = [profile_probability(text[slice(i, i+k)], profile) for i in xrange(n-k+1)]
    # randomly choose a starting position with probability proportional to probs
    j = random_choice(probs)
    # return chosen k-mer
    return text[slice(j, j+k)]

# one step of the gibbs sampler
# input:
#  k: k-mer length
#  t: number of strings in dna
#  N: number of gibbs sampler moves to try
#  dna: list of strings
# output:
#  a set lof k-mers from dna
def find_one_motif(k, t, N, dna):
    best_motifs = []
    motifs = []
    dna_len = len(dna[0])

    # generate a set of k-mers randomly
    for i in xrange(t):
        j = random.randrange(dna_len-k+1)
        kmer = dna[i][slice(j, j+k)]
        best_motifs.append(kmer)
        motifs.append(kmer)
    # score them
    best_score = motif_score(best_motifs)

    # try N steps
    for j in xrange(N):
        # choose a string in dna to resample
        i = random.randrange(t)

        # make profile from current k-mers except chosen string
        cur_motifs = motifs[:i] + motifs[(i+1):]
        old_motif = motifs[i]        
        profile = make_profile(cur_motifs, k)
        
        # now sample a new string with probability given by profile
        # sand score it
        motifs[i] = profile_sample(dna[i], profile, k)
        cur_score = motif_score(motifs)

        # keep this new string if it improves the score
        if cur_score < best_score:
            best_motifs = motifs
            best_score = cur_score
        else:
            # it didn't improve score, so switch back to original string
            motifs[i] = old_motif
    return best_motifs

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

    # run gibbs sampler once, and use the result as best motif so far
    best_motifs = find_one_motif(k, t, N, dna)
    best_score = motif_score(best_motifs)

    # run gibbs sampler multiple times
    for i in xrange(1, num_starts):
        cur_motifs = find_one_motif(k, t, N, dna)
        cur_score = motif_score(cur_motifs)

        # keep score for new motifs if they improve the score
        if cur_score < best_score:
            best_motifs = cur_motifs
            best_score = cur_score
    return best_motifs

filename = sys.argv[1]
with open(filename, 'r') as f:
    k, t, N = map(int, f.readline().strip().split())
    dna = []
    for i in xrange(t):
        dna.append(f.readline().strip())
    print '\n'.join(gibbs_sampling_search(k, t, N, dna))
