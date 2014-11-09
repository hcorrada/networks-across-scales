import sys
import numpy as np
import heapq

from Bio.SubsMat.MatrixInfo import pam250 as p250

# wrapper to get score from
# hash representing upper triangular matrix
def get_score(x, y, mat):
  return mat[(x,y)] if (x,y) in mat else mat[(y,x)]

# outputs the lcs alignment given
# backtrack 
def outputalign(backtrack, v, w, i, j, string1, string2):
    # uses this construct to avoid
    # deep recursion limit
    while True:
        if backtrack[i,j] == 'st' or (i==0 and j==0):
            return string1, string2
        if backtrack[i,j] == 'dn':
            i, string1, string2 = i-1, v[i-1] + string1, '-' + string2
        elif backtrack[i,j] == 'rt':
            j, string1, string2 = j-1, '-' + string1, w[j-1] + string2
        else:
            i, j, string1, string2 = i-1, j-1, v[i-1] + string1, w[j-1] + string2

# computes local alignment
# for strings v and w using dynamic programming
# and gap penalty sigma
def localalign(v, w, sigma):
    nv = len(v)
    nw = len(w)

    # initialize dp score and backtrack matrices
    s = np.zeros((nv+1,nw+1))
    backtrack = np.zeros((nv+1, nw+1), np.dtype('a2'))

    allscores = [(0, (0,0))]
    
    # fill in the dp matrices
    for i in xrange(0,nv+1):
        for j in xrange(0, nw+1):
            # use a heap to find max choice
            # use -score since heap keeps min in root

            # start with the 0 choice to stop the alignment
            choices = [(0, 'st')]

            if i > 0:
                heapq.heappush(choices, (-(s[i-1,j] + sigma), 'dn'))

            if j > 0:
                heapq.heappush(choices, (-(s[i,j-1] + sigma), 'rt'))

            if i>0 and j>0:
                score = get_score(v[i-1], w[j-1], p250)
                heapq.heappush(choices, (-(s[i-1,j-1] + score), 'di'))
            if len(choices) > 0:
                choice = choices[0]
                s[i,j] = -choice[0]
                heapq.heappush(allscores, (-s[i,j], (i,j)))
                backtrack[i,j] = choice[1]
    choice = allscores[0]
    score = -choice[0]
    i,j = choice[1]
    string1, string2 = outputalign(backtrack, v, w, i, j, '', '')
    return int(score), string1, string2

# read data
def readdat(filename):
    with open(filename, 'r') as f:
        v = f.readline().strip()
        w = f.readline().strip()
    return v,w

def main(filename):
    v,w  = readdat(filename)
    score, v,w = localalign(v, w, -5)
    print score
    print v
    print w

if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
