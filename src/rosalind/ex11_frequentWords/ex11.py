import sys

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
    # shorthand to get a k-mer
    ind = lambda x: slice(x,x+k)

    # see if this k-mer has been counted already
    checked = [False] * len(text)

    res = list()

    # current maximum count observed so far
    maxCount = 0

    # iterate over all k-mers
    for i in xrange(len(text)-k+1):
        # if this is counted already, skip it
        if not checked[i]:
            kmer = text[ind(i)]
            count = 1
            # check all k-mers that start in later positions
            for j in xrange((i+1), len(text)-k+1):
                # if it's not counted already and it matches the current
                # kmer, add it to the count (and mark it as checked)
                if not checked[j] and kmer == text[ind(j)]:
                    count = count + 1
                    checked[j] = True
            # if the new count is bigger than the current max
            # initialize the list of results
            if count > maxCount:
                res = [kmer]
                maxCount = count
            # if new count equals the current max, add it to list
            elif count == maxCount:
                res.append(kmer)
    return res

# get data filename as first argument to script
filename = sys.argv[1]
with open(filename, 'r') as f:
    text = f.readline().strip()
    k = int(f.readline().strip())
    kmers = mfk(text, k)
    print " ".join(kmers)
    
