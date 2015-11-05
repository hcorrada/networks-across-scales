from common import setup, readdat
from Bio.SubsMat.MatrixInfo import blosum62
import sys
import numpy as np

# find middle edge in a dynamic programming
# global alignment
# iv and iw are lists of indices into numpy scoring array
# see 'common.py' for details on how that works
#
# sigma and npmat are gap penalty and numpy scoring array
# top, bottom, left, right: define the portion of alignment being computed
#
# returns
#   source: coordinates of source of middle edge
#   sink: coordinates of sink of middle edge
#   score: score of maximal path going through middle edge
def middle_edge(iv, iw, sigma, npmat, top, bottom, left, right):
    middle = (left + right) / 2

    # the substring of v to align in this step (defined as a list of indices)
    # and its reverse
    cv = iv[top:bottom]
    rv = cv[::-1]

    # the first half of w (up to middle)
    # and the reverse of the second half of w
    cw = iw[left:middle]
    rw = iw[middle+1:right][::-1]

    # compute alignment scores from source (left) to middle column
    scoreleft, _ignored = scoreit(cv, cw, sigma, npmat)

    # compute alignment scores from sink (right) to middle column
    # this is done by reversing strings, so both results need to be
    # reversed afterwards
    # backtrack is a vector of arrows, indicating the kind of outgoing edge
    # for every node in the middle column
    scoreright, backtrack = scoreit(rv, rw, sigma, npmat)

    # reverse score from right since it was calculated on reversed strings
    scores = scoreleft + scoreright[::-1]

    # find the middle node and edge
    node = np.argmax(scores)
    edge = backtrack[::-1][node]
    print edge
    
    # and the score of the maximal path
    score = scores[node]

    # reconfigure where we are relative to complete graph
    top = top+node

    # this is the source node of the middle edge
    source = (top, middle)

    # get the sink node of the middle edge
    top += 0 if edge == 'rt' else 1
    middle += 0 if edge == 'dn' else 1
    sink = (top,middle)

    return source, sink, score

# score alignment between v and w (as lists of indices)
# gap penalty sigma and (numpy) scoring array npmat
# using linear space (workspace)
def scoreit(iv, iw, sigma, npmat):
    nv = len(iv)
    nw = len(iw)

    # allocate scoring workspace (only 2 columns)
    workspace = np.zeros((nv+1, 2))

    # allocate *vector* of backtrack arrows (we only care about the backtrack
    # arrows for the last column)
    backtrack = np.zeros((nv+1,1), np.dtype('a2'))

    # initialize first column of scoring workspace
    # using gap penalty
    workspace[:,0] = sigma * np.arange(0,nv+1)

    # now compute remaining columns
    for j in xrange(1, nw+1):
        # start with the right arrows
        workspace[:,1] = workspace[:,0] + sigma
        backtrack.fill('rt')

        # now the diagonal arrows
        scores = npmat[iv,iw[j-1]]

        # this is the score of moving along the diagonal for all nodes in the column
        tmp = workspace[:-1,0] + scores

        # check which of these improve on the right arrows
        ind = tmp > workspace[1:,1]

        # update entries in workspace and backtrack
        workspace[1:,1][ind] = tmp[ind]
        backtrack[1:][ind] = 'di'

        # now the down arrows
        for i in xrange(1,nv+1):
            tmp = workspace[i-1,1] + sigma
            if  tmp > workspace[i,1]:
                workspace[i,1] = tmp
                backtrack[i] = 'dn'

        if j < nw:
            # swap the two columns if this is not the last one
            workspace[:, 0] = workspace[:, 1]
    # return second column of workspace (since these were not swapped in the last step)
    # also return backtrack vector
    return workspace[:,1], backtrack

def main(filename):
    v, w = readdat(filename)
    iv, iw, npmat = setup(v, w, blosum62)
    source, sink, score = middle_edge(iv, iw, -5, npmat, 0, len(v), 0, len(w))
    print source, sink

if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
