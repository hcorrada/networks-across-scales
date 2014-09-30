import sys

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

def greedy_motif_search(k, t, dna):
    best_motifs = [string[:k] for string in dna]
    best_score = motif_score(best_motifs)
    
    first_string = dna[0]
    n = len(first_string)

    for i in xrange(n-k+1):
        motif_set = [first_string[slice(i, i+k)]]
        for j in xrange(1, t):
            profile = make_profile(motif_set, k)
            motif_set.append(most_probable_kmer(dna[j], k, profile))
        cur_score = motif_score(motif_set)
        if cur_score < best_score:
            best_score = cur_score
            best_motifs = motif_set
    return best_motifs
        

filename = sys.argv[1]
with open(filename, 'r') as f:
    k, t = map(int, f.readline().strip().split())
    dna = []
    for i in xrange(t):
        dna.append(f.readline().strip())
    print '\n'.join(greedy_motif_search(k, t, dna))
