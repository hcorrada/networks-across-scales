import sys
import numpy as np
import heapq

from Bio.SubsMat.MatrixInfo import blosum62 as bl62

# wrapper to get score from
# hash representing upper triangular matrix
def get_score(x, y, mat):
  return mat[(x,y)] if (x,y) in mat else mat[(y,x)]

# outputs the alignment given
# backtrack matrix, input strings and node in table
# to start from
def output_align(backtrack, v, w, i, j):
    v_alignment = ''
    w_alignment = ''
    
    # uses this construct to avoid
    # deep recursion limit
    while True:
        if i==0 and j==0:
            return v_alignment, w_alignment
        if backtrack[i,j] == 'dn':
            i, v_alignment, w_alignment = i-1, v[i-1] + v_alignment, '-' + w_alignment
        elif backtrack[i,j] == 'rt':
            j, v_alignment, w_alignment = j-1, '-' + v_alignment, w[j-1] + w_alignment
        else:
            i, j, v_alignment, w_alignment = i-1, j-1, v[i-1] + v_alignment, w[j-1] + w_alignment

# computes global alignment
# for strings v and w using dynamic programming
# and gap penalty sigma
def global_align(v, w, sigma):
    nv = len(v)
    nw = len(w)

    # initialize dp score and backtrack matrices
    s = np.zeros((nv+1,nw+1))
    backtrack = np.zeros((nv+1, nw+1), np.dtype('a2'))

    # fill in the dp matrices
    for i in xrange(0,nv+1):
        for j in xrange(0, nw+1):
            # use a heap to find max choice
            choices = []

            # use -score since heap keeps min in root
            if i > 0:
                heapq.heappush(choices, (-(s[i-1,j] + sigma), 'dn'))

            if j > 0:
                heapq.heappush(choices, (-(s[i,j-1] + sigma), 'rt'))

            if i>0 and j>0:
                score = get_score(v[i-1], w[j-1], bl62)
                heapq.heappush(choices, (-(s[i-1,j-1] + score), 'di'))

            if len(choices) > 0:
                choice = choices[0]
                s[i,j] = -choice[0]
                backtrack[i,j] = choice[1]

    score = s[nv, nw]
    return int(score), backtrack

# read data
def readdat(filename):
    with open(filename, 'r') as f:
        v = f.readline().strip()
        w = f.readline().strip()
    return v,w

def main(filename):
    v,w  = readdat(filename)
    score, backtrack = global_align(v, w, -5)
    v_alignment, w_alignment = output_align(backtrack, v, w, len(v), len(w))
    print score
    print v_alignment
    print w_alignment

if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
