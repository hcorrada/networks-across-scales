import numpy as np

# wrapper to get matrix score
def get_score(x, y, mat):
    return mat[(x,y)] if (x,y) in mat else mat[(y,x)]

def get_mat_alphabet(mat):
    return sorted(list(set(zip(*mat.keys())[0])))

def make_npmat(mat):
    alphabet = get_mat_alphabet(mat)
    n = len(alphabet)

    npmat = np.zeros((n,n))
    for i in xrange(n):
        for j in xrange(i,n):
            score = get_score(alphabet[i], alphabet[j], mat)
            npmat[i,j] = score
            npmat[j,i] = score
    return npmat

# read data from given file
def readdat(filename):
    with open(filename, 'r') as f:
        v = f.readline().strip()
        w = f.readline().strip()
    return v, w
    
def setup(v, w, mat):
    alphabet = get_mat_alphabet(mat)
    npmat = make_npmat(mat)

    iv = np.array(map(alphabet.index, list(v)))
    iw = np.array(map(alphabet.index, list(w)))
    return iv, iw, npmat

def check_result(v,w,va,wa,score,sigma,mat):
    assert len(va) == len(wa)
    assert v == va.replace('-','')
    assert w == wa.replace('-','')

    score2 = 0
    for i in xrange(len(va)):
        score2 += get_score(va[i],wa[i],mat) if va[i] != '-' and wa[i] != '-' else sigma
    assert score == score2
