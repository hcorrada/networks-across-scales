import sys
import numpy as np
import heapq

from Bio.SubsMat.MatrixInfo import pam250 as p250
v
# wrapper to get score from
# hash representing upper triangular matrix
def get_score(x, y, mat):
  return mat[(x,y)] if (x,y) in mat else mat[(y,x)]

# outputs the local alignment given
# backtrack starting at node (i,j)
def outputalign(backtrack, v, w, i, j):
    string1 = ''
    string2 =''
    
    # uses this construct to avoid
    # deep recursion limit
    while True:
        if backtrack[i,j] == 'st' or (i==0 and j==0):
            # starting over or reached source
            return string1, string2
        if backtrack[i,j] == 'dn':
            # took a down edge: gap in w
            i, string1, string2 = i-1, v[i-1] + string1, '-' + string2
        elif backtrack[i,j] == 'rt':
            # took a right edge: gap in v
            j, string1, string2 = j-1, '-' + string1, w[j-1] + string2
        else:
            # aligned, use a character from each string
            i, j, string1, string2 = i-1, j-1, v[i-1] + string1, w[j-1] + string2

# computes local alignment
# for strings v and w using dynamic programming
# and gap penalty sigma (usually a negative number)
def localalign(v, w, sigma):
    nv = len(v)
    nw = len(w)

    # initialize dp score and backtrack matrices
    s = np.zeros((nv+1,nw+1))

    # backtrack will contain codes for the type of arrow:
    # 'st': starting over (took the incoming edge from source node)
    # 'dn': took a "down" arrow: a gap in w
    # 'rt': took a "right" arrow: a gap in v
    # 'di': took a "diagonal" arrow: align character from v and w
    backtrack = np.zeros((nv+1, nw+1), np.dtype('a2'))

    best_score = -float('inf')
    best_node = (0,0)
    
    # fill in the dp matrices
    for i in xrange(0,nv+1):
        for j in xrange(0, nw+1):
            # use a heap to find max choice
            # use -score since heap keeps min in root

            # start with the 0 (incoming arrow from source node)
            choices = [(0, 'st')]

            if i > 0:
                # gap in w
                heapq.heappush(choices, (-(s[i-1,j] + sigma), 'dn'))

            if j > 0:
                # gap in v
                heapq.heappush(choices, (-(s[i,j-1] + sigma), 'rt'))

            if i>0 and j>0:
                # align characters
                score = get_score(v[i-1], w[j-1], p250)
                heapq.heappush(choices, (-(s[i-1,j-1] + score), 'di'))
                
            if len(choices) > 0:
                choice = choices[0]
                s[i,j] = -choice[0]
                backtrack[i,j] = choice[1]

                # update the best overall score, since we can start backtracking anywhere in graph
                if s[i,j] > best_score:
                    best_score = s[i,j]
                    best_node = (i,j)

    # start backtracking in best overall scoring node
    score = best_score
    i,j = best_node
    string1, string2 = outputalign(backtrack, v, w, i, j)
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
