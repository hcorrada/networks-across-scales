import sys

complement = dict(A="T", C="G", G="C", T="A")

def revcomp(text):
    out = list(text)
    out.reverse()
    for i in xrange(len(out)):
        out[i] = complement[out[i]]
    return ''.join(out)
    
filename = sys.argv[1]
with open(filename, 'r') as f:
    text = f.read().strip()
    res = revcomp(text)
    print res
