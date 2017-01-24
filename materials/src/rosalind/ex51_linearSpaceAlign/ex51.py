import sys
from common import setup, readdat, check_result
from middle_edge import middle_edge
import numpy as np
from Bio.SubsMat.MatrixInfo import blosum62

# compute global alignment for strings v and w
# using linear space dynamic programming
# with scoring matrix mat (a dictionary like those in Biopython) 
# and gap penalty sigma (usually a negative number)
def linearalign(v, w, mat, sigma):
    nv = len(v)
    nw = len(w)

    # to make things faster, we adjust the
    # scoring matrix to make it a symmetric numpy
    # array, we also convert strings v and w into list of indices
    # into the resulting numpy scoring array
    # see file 'common.py' for details
    iv, iw, npmat = setup(v, w, mat)
    
    # call the recursive solver aligning full pair of strings (index lists)
    # path will be a list of tuples, giving node coordinates for the maximal path
    path, score = helper(iv,iw, sigma, npmat, 0, nv, 0, nw)

    # construct alignment result
    # adding gaps where indicated by resulting path
    va, wa = makealign(v, w, path)

    # this function asserts a couple of properties
    # of the solution, defined in 'common.py'
    check_result(v,w,va,wa,score,sigma,mat)
    
    return score, va, wa

# this does the real work, iv and iw are lists of indices into the numpy scoring array
# sigma and npmat are gap penalty and numpy scoring array
# top, bottom: the range of rows to fill in
# left, right: the range of columns to fill in
def helper(iv, iw, sigma, npmat, top, bottom, left, right):
    # we catch a few base cases separately when
    # there is no 'middle' node to find

    # this checks if we are aligning an empty string
    if top == bottom or left == right:
        path, score = align_empty_string(iv, iw, sigma, npmat, top, bottom, left, right)
        return path, score

    # this checks if we are aligning a single character 
    if bottom-top == 1 or right-left == 1:
        path, score = align_one_character(iv, iw, sigma, npmat, top, bottom, left, right)
        return path, score

    # find the source and sink of the middle edge
    # along with the score of maximal path going through this edge
    source, sink, score = middle_edge(iv, iw, sigma, npmat, top, bottom, left, right)

    # divide and conquer!
    path1, _ignored = helper(iv, iw, sigma, npmat, top, source[0], left, source[1])
    path2, _ignored = helper(iv, iw, sigma, npmat, sink[0], bottom, sink[1], right)

    # put together resulting paths and return
    return path1 + path2, score

# base case when aligning an empty string
# this is called when either top==bottom or left==right
# here paths correspond to all right (or down) arrows
# score is given by aligning the required number of gaps
def align_empty_string(iv, iw, sigma, npmat, top, bottom, left, right):
    # check which kind this is
    # start assuming top==bottom
    kind = 'fat'
    if left == right:
        kind = 'tall'

    # how many gaps will be needed
    ngaps = (bottom-top) if kind == 'tall' else (right-left)
    score = sigma * ngaps

    # make the path along sideways (or down)
    if kind == 'fat':
        path = [(top,x) for x in xrange(left,right+1)]
    else:
        path = [(x,left) for x in xrange(top,bottom+1)]

    
    return path, score

# base case when aligning a single character
# this is called when either top + 1 == bottom or left +1 == right
# here paths correspond to all right (or down) arrows except for one diagonal or down arrow
# e.g.
# -> -> -> \ -> ->
# or
# -> -> -> -> -> -> |
#
# and similarily for down arrows
def align_one_character(iv, iw, sigma, npmat, top, bottom, left, right):
    # check which kind this is
    # start assuming top + 1 == bottom
    kind = 'fat'
    if right - left == 1:
        kind = 'tall'

    # score the alignment of the single character against each character on the other string
    scores = npmat[iv[top], iw[left:right]] if kind == 'fat' else npmat[iv[top:bottom], iw[left]]

    # and find the position that maximizes the alignment
    imax = np.argmax(scores)

    # how many gaps do we need to cover rest of the string
    ngaps = (right-left) if kind == 'fat' else (bottom-top)

    # now, check if it's better to align the single character or have two gaps instead
    if scores[imax] > 2*sigma:
        # better to align, so we need one fewer gap
        score = scores[imax] + sigma*(ngaps-1)

        # construct the path with diagonal edge at corresponding position
        if kind == 'fat':
            path = [(top,left+x) for x in xrange(imax+1)] + [(bottom,left+x) for x in xrange(imax+1,ngaps+1)]
        else:
            path = [(top+x,left) for x in xrange(imax+1)] + [(top+x,right) for x in xrange(imax+1,ngaps+1)]
    else:
        # we'll have two gaps instead of an alignment
        # so, we'll need an extra gap
        score = sigma*(ngaps+1)
        if kind == 'fat':
            # walk along to the right, and add one down edge
            path = [(top,x) for x in xrange(left,right+1)] + [(bottom,right)]
        else:
            # walk along to the bottom, and add one right edge
            path = [(x,left) for x in xrange(top,bottom+1)] + [(bottom,right)]
    return path, score
        
# construct alignment from path
# path is list of pairs with coordinates
# in the graph where the path goes through
def makealign(v, w, path):
    string1 = ''
    string2 = ''
    i,j = path[0]
    for k in xrange(1, len(path)):
        # look at the next node
        i2, j2 = path[k]

        # gap or alignment: check if row coordinate changes
        string1 += v[i2-1] if i2 != i else '-'

        # gap or alignment: check if column coordinate changes
        string2 += w[j2-1] if j2 != j else '-'

        # update current indices        
        i, j = i2, j2
    return string1, string2

def main(filename):
    v, w = readdat(filename)
    score, v, w = linearalign(v, w, blosum62, -5)
    print int(score)
    print v
    print w
    
if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
