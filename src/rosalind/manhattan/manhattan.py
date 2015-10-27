import sys
import numpy as np

def manhattan(n, m, down_mat, right_mat):
    return None

def readdat(filename):
    with open(filename, 'r') as f:
        n, m = [int(s.strip()) for s in f.readline().split()]
        down_mat = np.zeros((n, m+1))
        right_mat = np.zeros((n+1, m))

        for i in xrange(n):
            a = np.array([int(x) for x in f.readline().split()])
            down_mat[i,] = a

        f.readline()

        for i in xrange(n+1):
            a = np.array([int(x) for x in f.readline().split()])
            right_mat[i,] = a

        return n, m, down_mat, right_mat

def main(filename):
    n, m, down_mat, right_mat = readdat(filename)
    print n, m
    print down_mat
    print right_mat
    
    path_length = manhattan(n, m, down_mat, right_mat)
    print path_length

if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
