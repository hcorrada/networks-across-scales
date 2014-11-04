import sys
import numpy as np
import heapq

def outputlcs(backtrack, v, i, j, lcs):
    while True:
        if i==0 or j==0:
            return lcs
        if backtrack[i,j] == 'dn':
            i = i-1
        elif backtrack[i,j] == 'rt':
            j = j-1
        else:
            i, j, lcs = i-1, j-1, v[i-1] + lcs

def lcs(v,w):
    nv = len(v)
    nw = len(w)
    
    s = np.zeros((nv+1,nw+1))
    backtrack = np.zeros((nv+1,nw+1), np.dtype('a2'))

    for i in xrange(0,nv+1):
        s[i,0] = 0
    for j in xrange(0, nw+1):
        s[0,j] = 0

    for i in xrange(1,nv+1):
        for j in xrange(1, nw+1):
            choices = []
            heapq.heappush(choices, (-s[i-1,j], 'dn'))
            heapq.heappush(choices, (-s[i,j-1], 'rt'))
            if (v[i-1] == w[j-1]):
                heapq.heappush(choices, (-(s[i-1,j-1] + 1), 'di'))
            choice = choices[0]
            s[i,j] = -choice[0]
            backtrack[i,j] = choice[1]
    out = outputlcs(backtrack, v, nv, nw, '')
    return out

def readdat(filename):
    with open(filename, 'r') as f:
        v = f.readline().strip()
        w = f.readline().strip()
    return v,w

def main(filename):
    v,w  = readdat(filename)
    print lcs(v, w)

if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
