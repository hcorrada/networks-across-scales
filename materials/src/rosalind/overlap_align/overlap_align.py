import sys
import numpy as np
import heapq

# wrapper to get score from
# hash representing upper triangular matrix
def get_score(x, y, mat):
    out = 1 if x==y else -2
    return out


# outputs the local alignment given
# backtrack starting at node (i,j)
def output_align(backtrack, v, w, i, j):
    string1 = ''
    string2 = ''

    # uses this construct to avoid
    # deep recursion limit
    while True:
        # I used this to debug
        # print string1, string2, i, j
        # HERE WAS THE TYPO (WE WROTE i==0 IN CLASS)
        if j==0:
            # reached first column
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

# computes global alignment
# for strings v and w using dynamic programming
# and gap penalty sigma
def overlap_align(v, w):
    sigma = -2

    nv = len(v)
    nw = len(w)

    # initialize dp score and backtrack matrices
    s = np.zeros((nv+1,nw+1))
    backtrack = np.zeros((nv+1, nw+1), np.dtype('a2'))

    # initialize dp table
    for i in xrange(0, nv+1):
        s[i,0] = 0
    for j in xrange(0, nw+1):
        s[0,j] = j * sigma

    best_score = s[nv,0]
    best_node = (nv, 0)

    # fill in the dp matrices
    for j in xrange(1, nw+1):
        for i in xrange(1, nv+1):
            # use a heap to find max choice
            choices = []

            # use -score since heap keeps min in root
            heapq.heappush(choices, (-(s[i-1, j] + sigma), 'dn'))

            heapq.heappush(choices, (-(s[i, j-1] + sigma), 'rt'))

            score = get_score(v[i-1], w[j-1], None)
            heapq.heappush(choices, (-(s[i-1,j-1] + score), 'di'))

            if len(choices) > 0:
                choice = choices[0]
                s[i,j] = -choice[0]
                backtrack[i,j] = choice[1]

        if s[nv, j] > best_score:
            best_score = s[nv, j]
            best_node = (nv, j)

    return int(best_score), backtrack, best_node

# read data
def readdat(filename):
    with open(filename, 'r') as f:
        v = f.readline().strip()
        w = f.readline().strip()
    return v,w

def check_solution(v,va,w,wa,score):
    assert len(va) == len(wa)
    score2 = 0
    for i in xrange(len(va)):
        score2 += 1 if va[i] == wa[i] else -2
    assert score2 == score

def main(filename):
    v,w  = readdat(filename)
    score, backtrack, (i,j) = overlap_align(v, w)
    va, wa = output_align(backtrack, v, w, i, j)
    check_solution(v,va,w,wa,score)
    print score
    print va
    print wa

if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
