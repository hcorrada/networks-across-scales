import sys
import random

random.seed(2)

alphabet = list('ACGT')

def motif_score(motifs):
    k = len(motifs[0])
    t = len(motifs)
    
    score = 0
    for i in xrange(k):
        count = dict.fromkeys(alphabet, 0)
        for j in xrange(t):
            count[motifs[j][i]] += 1
        most_common = max(count.items(), key = lambda x: x[1])[0]
        cur_score = sum(count.values()) - count[most_common]
        score += cur_score
    return score

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

def random_choice(weights):
    total_weight = sum(weights)
    nitems = len(weights)
    for i in xrange(nitems):
        weights[i] /= total_weight
        
    choice = -1
    u = random.uniform(0, 1.)
    
    while u >= 0:
        choice += 1
        u -= weights[choice]
    return choice

def profile_probability(kmer, profile):
    prob = 1.
    k = len(kmer)
    for i in xrange(k):
        prob *= profile[kmer[i]][i]
    return prob

def profile_sample(text, profile, k):
    n = len(text)
    probs = [profile_probability(text[slice(i, i+k)], profile) for i in xrange(n-k+1)]
    j = random_choice(probs)
    return text[slice(j, j+k)]

def find_one_motif(k, t, N, dna):
    best_motifs = []
    motifs = []
    dna_len = len(dna[0])
    
    for i in xrange(t):
        j = random.randrange(dna_len-k+1)
        kmer = dna[i][slice(j, j+k)]
        best_motifs.append(kmer)
        motifs.append(kmer)
    best_score = motif_score(best_motifs)

    for j in xrange(N):
        i = random.randrange(t)
        cur_motifs = motifs[:i] + motifs[(i+1):]
        old_motif = motifs[i]
        
        profile = make_profile(cur_motifs, k)
        motifs[i] = profile_sample(dna[i], profile, k)
        cur_score = motif_score(motifs)
        if cur_score < best_score:
            best_motifs = motifs
            best_score = cur_score
        else:
            motifs[i] = old_motif
    return best_motifs

def gibbs_sampling_search(k, t, N, dna):
    num_starts = 20
    best_motifs = find_one_motif(k, t, N, dna)
    best_score = motif_score(best_motifs)

    for i in xrange(1, num_starts):
        cur_motifs = find_one_motif(k, t, N, dna)
        cur_score = motif_score(cur_motifs)

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
