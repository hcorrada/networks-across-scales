import sys
from collections import defaultdict

# use this hash table to get complementary nucleotides
complement = dict(A="T", C="G", G="C", T="A")

# compute the reverse complement of a given string
# input:
#  text: a string
# output:
#  the reverse complement of text
def revcomp(text):
    # turn the string into a list
    out = list(text)

    # reverse the list (in place)
    out.reverse()

    # get complement at each position of list
    for i in xrange(len(out)):
        out[i] = complement[out[i]]

    # turn back into string and return
    return ''.join(out)

# generate the d neighborhood of a DNA string
nucleotides = list('ACGT')

def get_neighbors(pattern, d):
    if d == 0:
        return [(pattern,0)]
    if len(pattern) == 1:
        return [(nuc, 1 if nuc != pattern else 0) for nuc in nucleotides]

    neighborhood = set()
    suffix_neighbors = get_neighbors(pattern[1:], d)
    for (text, cur_d) in suffix_neighbors:
        if cur_d < d:
            for nuc in nucleotides:
                distance = cur_d + 1 if nuc != pattern[0] else cur_d
                neighborhood.add((nuc + text, distance))
        else:
            neighborhood.add((pattern[0] + text, d))
    return neighborhood

# compute most frequent words allowing d mismatches, including reverse complements
# input:
#   text: input string
#   k: k-mer length
#   d: maximum number of mismatches
def freq_words_mismatch_revcomp(text, k, d):
    kmer_counts = defaultdict(int)

    n = len(text)
    for i in xrange(n-k+1):
        kmer = text[slice(i,i+k)]
        neighbors = get_neighbors(kmer, d)
        for neighbor, _ in neighbors:
            kmer_counts[neighbor] += 1
        kmer = revcomp(kmer)
        neighbors = get_neighbors(kmer, d)
        for neighbor, _ in neighbors:
            kmer_counts[neighbor] += 1

    max_count = 0
    frequent_words = list()
    for kmer, count in kmer_counts.items():
        if count == max_count:
            frequent_words.append(kmer)
        if count > max_count:
            frequent_words = [kmer]
            max_count = count
    return frequent_words, max_count

filename = sys.argv[1]
with open(filename, 'r') as f:
    text = f.readline().strip()
    k, d = map(int, f.readline().strip().split())
    frequent_words, _ = freq_words_mismatch_revcomp(text, k, d)
    print " ".join(frequent_words)
