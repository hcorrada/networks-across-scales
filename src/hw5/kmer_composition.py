import sys

def kmer_composition(k, text):
    kmers = []
    for i in xrange(len(text) - k +1):
        kmers.append(text[slice(i,i+k)])
    kmers.sort()
    return kmers

def readdat(filename):
    with open(filename, 'r') as f:
        k = int(f.readline().strip())
        test = f.readline().strip()
    return k, test

def main(filename):
    k, text = readdat(filename)
    kmers = kmer_composition(k, text)
    print "\n".join(kmers)
    
if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
    
    
