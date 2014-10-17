import itertools

complement = dict(A="T", C="G", G="C", T="A")

def revcomp(text):
    out = list(text)
    out.reverse()
    for i in xrange(len(out)):
        out[i] = complement[out[i]]
    return ''.join(out)

def compute(n,k):
    total = math.pow(4, n)
    allseqs = itertools.product('ACGT', repeat=n)
    ok = 0
    
    for seq in allseqs:
        x = ''.join(seq)
        if revcomp(x[:k]) == x[-k:]:
            ok += 1
    return float(ok) / total

n = 8
res = [compute(n, k) for k in xrange(1,n)]
