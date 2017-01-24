import sys

# count the number of times pattern occurs in text
# input:
#   pattern: a string
#   text: a string
def count_pattern(pattern, text):
    k = len(pattern)
    n = len(text)

    count = 0
    ind = lambda x: slice(x, x+k)
    for i in xrange(n-k+1):
        if pattern == text[ind(i)]:
            count += 1
    return count

# find the most frequent k-mers in a given string
# input:
#  text: string
#  k: k-mer length
#
# output:
#  list of most frequent k-mers
#
# this is a quadratic algorithm
def mfk(text, k):
    n = len(text)

    # count the number of times kmer at position i
    # occurs in text
    counts = [0] * (n-k+1)
    for i in xrange(len(counts)):
        kmer = text[slice(i, i+k)]
        counts[i] = count_pattern(kmer, text)

    # now find the maximum kmers
    # current maximum count observed so far
    maxCount = 0

    # iterate over all k-mers
    for i in xrange(len(counts)):
            kmer = text[slice(i,i+k)]
            if counts[i] > maxCount:
                res = [kmer]
                maxCount = counts[i]
            # if new count equals the current max, add it to list
            elif counts[i] == maxCount and kmer not in res:
                res.append(kmer)
    return res

# get data filename as first argument to script
filename = sys.argv[1]
with open(filename, 'r') as f:
    text = f.readline().strip()
    k = int(f.readline().strip())
    kmers = mfk(text, k)
    print " ".join(kmers)
