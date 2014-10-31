import sys
import random

def kdcomposition(text, k, d):
    total_len = 2*k+d
    n = len(text)
    read_pairs = []
    for i in xrange(n-total_len+1):
        kmer1=text[slice(i,i+k)]
        kmer2=text[slice(i+k+d,i+total_len)]
        read_pairs.append(kmer1 + "|" + kmer2)
    random.shuffle(read_pairs)
    return read_pairs

k = 50
d = 200

filename = sys.argv[1]
with open(filename,'r') as f:
    text = f.read().strip()
    read_pairs = kdcomposition(text, 50, 200)
    print d
    print "\n".join(read_pairs)
