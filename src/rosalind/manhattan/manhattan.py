import sys
import numpy as np

def manhattan(n, m, down_mat, right_mat):
    dp_table = np.zeros((n+1, m+1))

    # intialize left column
    for i in xrange(1, n+1):
        dp_table[i,0] = dp_table[i-1, 0] + down_mat[i-1, 0]

    # initialize top row
    for j in xrange(1, m+1):
        dp_table[0,j] = dp_table[0,j-1] + right_mat[0, j-1]

    # fill-in the rest of the table, row by row
    for i in xrange(1, n+1):
        for j in xrange(1, m+1):
            val_down = dp_table[i-1,j] + down_mat[i-1,j]
            val_right = dp_table[i, j-1] + right_mat[i,j-1]
            choices = np.array([val_down, val_right])
            choice = np.argmax(choices)
            dp_table[i,j] = choices[choice]

    return dp_table[n, m]

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
    path_length = manhattan(n, m, down_mat, right_mat)
    print int(path_length)

if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
