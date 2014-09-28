import sys
from collections import deque

def clump_finder(genome, k, L, t):
    ind = lambda x: slice(x, x+k)
    checked = [False] * len(genome)

    res = list()
    
    for i in xrange(len(genome) - k + 1):
        if not checked[i]:
            kmer = genome[ind(i)]
            occurences = deque([i])
            done = False

            for j in xrange((i+1), len(genome) - k +1):
                if not checked[j] and kmer == genome[ind(j)]:
                    checked[j] = True
                    if done:
                        next
                    while len(occurences) > 0 and occurences[0] <= j - L + 1:
                        occurences.popleft()
                    occurences.append(j)
                    
            if len(occurences) >= t:
                res.append(kmer)
                done = True
    return res

filename = sys.argv[1]
with open(filename, 'r') as f:
    genome = f.readline().strip()
    k, L, t = map(int, f.readline().strip().split())
    res = clump_finder(genome, k, L, t)
    print ' '.join(res)
