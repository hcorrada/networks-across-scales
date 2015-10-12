import sys

# compute the k-mer composition of string
# input:
#   k: k-mer size
#   text: dna string
# output:
#   list containing all k-mers in text
def kmer_composition(k, text):
    n = len(text)
    kmers = []
    for i in xrange(n - k +1):
        kmers.append(text[slice(i,i+k)])
    # sort in place
    kmers.sort()
    return kmers

# read input from file
# expected format:
#   k
#   text
# output:
#   tuple k, text
def readdat(filename):
    with open(filename, 'r') as f:
        k = int(f.readline().strip())
        test = f.readline().strip()
    return k, test

# read data from file, print the k-mer composition
def main(filename):
    k, text = readdat(filename)
    kmers = kmer_composition(k, text)
    print "\n".join(kmers)

# this is here to make it easier to work with ipython
if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
