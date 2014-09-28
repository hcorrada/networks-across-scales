import sys

def mfk(text, k):
    ind = lambda x: slice(x,x+k)
    checked = [False] * len(text)

    res = list()
    maxCount = 0
    
    for i in xrange(len(text)-k+1):
        if not checked[i]:
            kmer = text[ind(i)]
            count = 1
            for j in xrange((i+1), len(text)-k+1):
                if not checked[j] and kmer == text[ind(j)]:
                    count = count + 1
                    checked[j] = True
            if count > maxCount:
                res = [kmer]
                maxCount = count
            elif count == maxCount:
                res.append(kmer)
    return res

filename = sys.argv[1]
with open(filename, 'r') as f:
    text = f.readline().strip()
    k = int(f.readline().strip())
    kmers = mfk(text, k)
    print " ".join(kmers)
    
