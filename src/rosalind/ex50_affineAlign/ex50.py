import sys
import numpy as np
import heapq

from Bio.SubsMat.MatrixInfo import blosum62 as bl62

# wrapper to get score from
# hash representing upper triangular matrix
def get_score(x, y, mat):
  return mat[(x,y)] if (x,y) in mat else mat[(y,x)]

# outputs the lcs alignment given
# backtrack 
def outputalign(backtrack_lower, backtrack_middle, backtrack_upper, v, w, curmat, i, j, string1, string2):
    # uses this construct to avoid
    # deep recursion limit
    while True:
        if i==0 and j==0:
            return string1, string2
        
        if curmat == 'lo':
            i, string1, string2, curmat = i-1, v[i-1] + string1, '-' + string2, backtrack_lower[i,j]
        elif curmat == 'up':
            j, string1, string2, curmat = j-1, '-' + string1, w[j-1] + string2, backtrack_upper[i,j]
        elif backtrack_middle[i,j] != 'mi':
            curmat = backtrack_middle[i,j]
        else:
            i, j, string1, string2, curmat = i-1, j-1, v[i-1] + string1, w[j-1] + string2, backtrack_middle[i,j]

# computes global affine alignment
# for strings v and w using dynamic programming
# and gap open penalty sigma and gap extension epsilon
def affinealign(v, w, sigma, epsilon):
    nv = len(v)
    nw = len(w)

    # initialize dp score and backtrack matrices
    s_lower = np.zeros((nv+1, nw+1))
    s_middle = np.zeros((nv+1, nw+1))
    s_upper = np.zeros((nv+1, nw+1))
    
    backtrack_lower = np.zeros((nv+1, nw+1), np.dtype('a2'))
    backtrack_middle = np.zeros((nv+1, nw+1), np.dtype('a2'))
    backtrack_upper = np.zeros((nv+1, nw+1), np.dtype('a2'))

    # fill in the dp matrices
    for i in xrange(0,nv+1):
        for j in xrange(0, nw+1):

            # use a heap to find max choice
            # use -score since heap keeps min in root
            if i > 0:
                choices = []
                heapq.heappush(choices, (-(s_lower[i-1,j] + epsilon), 'lo'))
                heapq.heappush(choices, (-(s_middle[i-1,j] + sigma), 'mi'))
                choice = choices[0]
                s_lower[i,j] = -choice[0]
                backtrack_lower[i,j] = choice[1]                
            if j > 0:
                choices = []
                heapq.heappush(choices, (-(s_upper[i,j-1] + epsilon), 'up'))
                heapq.heappush(choices, (-(s_middle[i,j-1] + sigma), 'mi'))
                choice = choices[0]
                s_upper[i,j] = -choice[0]
                backtrack_upper[i,j] = choice[1]
            if i>0 and j>0:
                choices = []
                score = get_score(v[i-1], w[j-1], bl62)
                heapq.heappush(choices, (-(s_lower[i,j]), 'lo'))
                heapq.heappush(choices, (-(s_middle[i-1,j-1] + score), 'mi'))
                heapq.heappush(choices, (-(s_upper[i,j]), 'up'))
                choice = choices[0]
                s_middle[i,j] = -choice[0]
                backtrack_middle[i,j] = choice[1]

    # choose which of the three matrices to start the alignment from
    choices = []
    heapq.heappush(choices, (-s_lower[nv,nw], 'lo'))
    heapq.heappush(choices, (-s_middle[nv,nw], 'mi'))
    heapq.heappush(choices, (-s_upper[nv,nw], 'up'))
    choice = choices[0]
    score = -choice[0]
    curmat = choice[1]

    string1, string2 = outputalign(backtrack_lower, backtrack_middle, backtrack_upper, v, w, curmat, nv, nw, '', '')
    return int(score), string1, string2

# read data
def readdat(filename):
    with open(filename, 'r') as f:
        v = f.readline().strip()
        w = f.readline().strip()
    return v,w

def main(filename):
    v,w  = readdat(filename)
    score, v,w = affinealign(v, w, -11, -1)
    print score
    print v
    print w

if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
