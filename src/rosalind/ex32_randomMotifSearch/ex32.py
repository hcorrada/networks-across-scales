import sys
import random

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

def profile_probability(kmer, profile):
    prob = 1.
    k = len(kmer)
    for i in xrange(k):
        prob *= profile[kmer[i]][i]
    return prob

def most_probable_kmer(text, k, profile):
    best_kmer = text[:k]
    best_probability = profile_probability(best_kmer, profile)

    n = len(text)
    for i in xrange(1, n-k+1):
        kmer = text[slice(i,i+k)]
        prob = profile_probability(kmer, profile)
        if prob > best_probability:
            best_probability = prob
            best_kmer = kmer
    return best_kmer

def make_motifs(profile, dna, t, k):
    motifs = []
    for i in xrange(t):
        motifs.append(most_probable_kmer(dna[i], k, profile))
    return motifs

def find_one_motif(k, t, dna):
    best_motifs = []
    motifs = []
    n = len(dna[0])
    
    for i in xrange(t):
        j = random.randrange(n-k+1)
        kmer = dna[i][slice(j, j+k)]
        best_motifs.append(kmer)
        motifs.append(kmer)
    best_score = motif_score(best_motifs)

    while True:
        profile = make_profile(motifs, k)
        motifs = make_motifs(profile, dna, t, k)
        cur_score = motif_score(motifs)
        if cur_score < best_score:
            best_motifs = motifs
            best_score = cur_score
        else:
            return best_motifs

def random_motif_search(k, t, dna):
    n = 1000
    best_motifs = find_one_motif(k, t, dna)
    best_score = motif_score(best_motifs)
    for i in xrange(1, n):
        cur_motifs = find_one_motif(k, t, dna)
        cur_score = best_score
        if cur_score < best_score:
            best_motifs = cur_motifs
            cur_score = best_score
    return best_motifs

filename = sys.argv[1]
with open(filename, 'r') as f:
    k, t = map(int, f.readline().strip().split())
    dna = []
    for i in xrange(t):
        dna.append(f.readline().strip())
    print '\n'.join(random_motif_search(k, t, dna))
