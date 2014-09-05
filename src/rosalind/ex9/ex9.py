import sys

complement = dict(A="T", C="G", G="C", T="A")
        
filename = sys.argv[1]
with open(filename, 'r') as f:
    out = list(f.read().strip())
    out.reverse()
    for i in xrange(len(out)):
        out[i] = complement[out[i]]
    print ''.join(out)
