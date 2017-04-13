import sys
import numpy as np
import heapq

from Bio.SubsMat.MatrixInfo import blosum62 as bl62

# wrapper to get score from
# hash representing upper triangular matrix
def get_score(x, y, mat):
  return mat[(x,y)] if (x,y) in mat else mat[(y,x)]

# outputs the alignment given
# backtrack matrix, input strings and position in table
# to start backtrack from
def output_align(backtrack, a, b, i, j):
    a_alignment = ''
    b_alignment = ''
    
    # uses this construct to avoid
    # deep recursion limit
    while True:
        if i==0 and j==0:
            return a_alignment, b_alignment
        if backtrack[i,j] == 'dn':
            i, a_alignment, b_alignment = i-1, a[i-1] + a_alignment, '-' + b_alignment
        elif backtrack[i,j] == 'rt':
            j, a_alignment, b_alignment = j-1, '-' + a_alignment, b[j-1] + b_alignment
        else:
            i, j, a_alignment, b_alignment = i-1, j-1, a[i-1] + a_alignment, b[j-1] + b_alignment

# computes global alignment
# for strings v and w using dynamic programming
# and gap penalty sigma
def global_align(a, b, sigma):
    n = len(a)
    m = len(b)

    # initialize dp score and backtrack matrices
    s = np.zeros((n+1,m+1))
    backtrack = np.zeros((n+1, m+1), np.dtype('a2'))

    # initialize dp table
    for i in xrange(0, n+1):
        s[i,0] = -sigma * i
        backtrack[i,0] = 'dn'

    for j in xrange(0, m+1):
        s[0,j] = -sigma * j
        backtrack[0,j] = 'rt'

    # fill in the dp matrices
    for i in xrange(1, n+1):
        for j in xrange(1, m+1):
            # use a heap to find max choice
            choices = []

            # use -score since heap keeps min in root
            heapq.heappush(choices, (-(s[i-1,j] - sigma), 'dn'))
            heapq.heappush(choices, (-(s[i,j-1] - sigma), 'rt'))

            score = get_score(a[i-1], b[j-1], bl62)
            heapq.heappush(choices, (-(s[i-1,j-1] + score), 'di'))

            choice = choices[0]
            s[i,j] = -choice[0]
            backtrack[i,j] = choice[1]

    score = s[n, m]
    return int(score), backtrack

# read data
def readdat(filename):
    with open(filename, 'r') as f:
        a = f.readline().strip()
        b = f.readline().strip()
    return a, b

def main(filename):
  sigma = 5 # given in problem definition
  a, b  = readdat(filename)
  score, backtrack = global_align(a, b, sigma)
  a_alignment, b_alignment = output_align(backtrack, a, b, len(a), len(b))
  print score
  print a_alignment
  print b_alignment

if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
