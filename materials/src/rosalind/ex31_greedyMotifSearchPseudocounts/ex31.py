import sys

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

# compute the probabilty of a given k-mer according to given profile
# input:
#  kmer: a k-mer
#  profile: a profile
# output:
#  the probability
def profile_probability(kmer, profile):
    # initialize probability
    prob = 1.
    k = len(kmer)
    for i in xrange(k):
        # for each character in kmer, multiply in the probability of that nucleotide
        # according to profile
        prob *= profile[kmer[i]][i]
    return prob

# find the most probable k-mer in a given string according to given profile
# input:
#  text: a string
#  k: k-mer length
#  profile: a profile
# output
#   the most probable k-mer in text
def most_probable_kmer(text, k, profile):
    # start with the first k-mer in text
    best_kmer = text[:k]
    best_probability = profile_probability(best_kmer, profile)

    n = len(text)
    # compute probability of each k-mer in text
    for i in xrange(1, n-k+1):
        kmer = text[slice(i,i+k)]
        prob = profile_probability(kmer, profile)

        # keep this k-mer if most probable so far
        if prob > best_probability:
            best_probability = prob
            best_kmer = kmer
    return best_kmer

# greedy motif search
# input:
#  k: k-mer length
#  t: number of strings in dna
#  dna: a list of DNA strings
# output:
#  the best set of motifs found
def greedy_motif_search(k, t, dna):
    # initialize by selecting the first k-mer in each string
    best_motifs = [string[:k] for string in dna]
    best_score = motif_score(best_motifs)
    
    first_string = dna[0]
    n = len(first_string)

    # get the best solution greedily, selecting each k-mer in first string
    for i in xrange(n-k+1):
        # initialize this current motif set with k-mer in first string
        motif_set = [first_string[slice(i, i+k)]]

        # add k-mers from other strings in order
        for j in xrange(1, t):
            # make a profile with current set of strings
            profile = make_profile(motif_set, k)

            # add the most-probable k-mer in new string to current set
            motif_set.append(most_probable_kmer(dna[j], k, profile))
        # if the resulting set is the best so far, keep it
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
