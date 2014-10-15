import sys
import itertools

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

# find the kmers that occur in a string and compute the number of times they appear
# input:
#  text: a string
#  k: k-mer length
# output:
#   a hash table indexed by k-mer with the number of occurences as the entry
def count_kmers(text, k):
    n = len(text)

    # start and empty hash table
    tab = {}

    # for every position in text (where a k-mer can start)
    for i in xrange(n-k+1):
        # get the k-mer
        kmer = text[slice(i, i+k)]
        # add 1 to the current count if it's in the table, otherwise add to table with count of 1
        tab[kmer] = tab[kmer] + 1 if kmer in tab else 1
    return tab

# given a k-mer pattern, return the list of all k-mers within d mismatches
# input:
#   pattern: a k-mer
#   d: maximum number of mismatches
# output:
#   a list of k-mers within d mismatches of pattern
def build_neighborhood(pattern, d):
    # use a hash table to avoid adding a k-mer more than once
    res = {}
    k = range(len(pattern))

    # select d positions out of the k
    for i in itertools.combinations(k, d):
        # generate all possible d-combinations of nucleotides
        for j in itertools.product('ACGT', repeat=d):
            # replace the generated d-combination in the selected positions
            # turn the string into a list
            temp=list(pattern)
            for k in xrange(d):
                temp[i[k]] = j[k]
            # turn back into string
            temp = "".join(temp)
            # add to hash table if it's not there
            if temp not in res:
                res[temp] = 0
    return res.keys()

# process a k-mer (see function below for explanation)
# input:
#  pattern: a kmer
#  final_tab: a table containing running counts for each candidate k-mer
#  kmer_counts: the k-mer count table from input
#  d: maximum number of mismatches
def process_kmer(pattern, final_tab, kmer_counts, d):
    # generate all k-mers within d mismatches of patterns
    neighbors = build_neighborhood(pattern, d)

    # for each of these k-mers
    for kmer in neighbors:
        # if it's already in the table, add to it the number of times pattern appears in the input string
        if kmer in final_tab:
            final_tab[kmer] += kmer_counts[pattern]
        # otherwise add it to the final hash table
        else:
            final_tab[kmer] = kmer_counts[pattern]
    return final_tab

# compute most frequent words allowing d mismatches, including reverse complements
# input:
#   text: input string
#   k: k-mer length
#   d: maximum number of mismatches
def freq_words_mismatch(text, k, d):
    # count k-mers in input string
    kmer_counts = count_kmers(text, k)

    # get the list of kmers out of the hash table
    kmers = kmer_counts.keys()

    # initialize hash table
    final_tab = {}
    for kmer in kmers:
        # process each k-mer: generate all k-mers within d mismatches
        #   add the number of times kmer occurs in the input to the running count
        #   stored in final_tab
        final_tab = process_kmer(kmer, final_tab, kmer_counts, d)

    # now go through the table and find the maximum k-mers
    maxcount = 0
    res = []

    for (kmer, count) in final_tab.iteritems():
        rc_kmer = revcomp(kmer)
        if rc_kmer in final_tab:
            count += final_tab[rc_kmer]
        
        if count > maxcount:
            maxcount = count
            res = [kmer]
        elif count == maxcount:
            res.append(kmer)
    return res
    
filename = sys.argv[1]
with open(filename, 'r') as f:
    text = f.readline().strip()
    k, d = map(int, f.readline().strip().split())
    print " ".join(freq_words_mismatch(text, k, d))
    
