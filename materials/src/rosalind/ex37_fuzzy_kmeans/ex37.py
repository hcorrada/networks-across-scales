import sys
import heapq
import math

# compute euclidean distance between
# vectors x and y of length m
def point_to_point_distance(x, y, m):
    res = 0
    for i in xrange(m):
        res += math.pow(x[i] - y[i], 2)
    return math.sqrt(res)

# return the minimum euclidean distance
# between point y and set of points x
# points are of dimension m
def point_to_set_distance(y, x, m):
    heap = []
    for i in xrange(len(x)):
        heapq.heappush(heap, point_to_point_distance(y, x[i], m))
    return heap[0]

# initialize the set of centers
# using furthest traversal
def initialize(Data, k, n, m):
    # choose the first point in Data as a center
    centers = [0]

    # initialize the matrix of returned centers
    x = [Data[0]]

    # while more centers are needed
    while len(centers) < k:
        # use a heap to find next center:
        # remaining point in Data that maximizes the distance
        # to the current set of centers
        dist = []
        for i in xrange(n):
            # this point is already a center, skip it
            if i in centers:
                continue
            # push item to heap
            # item is (-distance to centers, row index in Data)
            # use -distance since heap keeps the smallest value at root
            heapq.heappush(dist, (-point_to_set_distance(Data[i], x, m), i))

        # next center is at root of the heap
        # item on heat is (-distance, index)
        next_center = dist[0][1]

        # add the center to set
        centers.append(next_center)
        x.append(Data[next_center])
    return x

# compute hidden responsibility vector for a point
def get_hidden_vector(y, X, beta, m):
    vec = []
    denom = 0.
    
    for i in xrange(len(X)):
        d = point_to_point_distance(y, X[i], m)
        d = math.exp(-beta * d)
        denom += d
        vec.append(d)
    for i in xrange(len(X)):
        vec[i] /= denom
    return vec

# compute hidden matrix of responsibilities
def get_hidden_matrix(Data, X, beta, n, m):
    k = len(X)
    hidden_matrix = []
    for i in xrange(n):
        hidden_matrix.append(get_hidden_vector(Data[i], X, beta, m))
    return hidden_matrix

# get k-th center of a set of points
# weighted by k-th column in hidden_matrix
def get_center(X, hidden_matrix, k, m):
    center = []
    n = len(X)
    
    for j in xrange(m):
        running_sum = 0.
        denom = 0.
        
        for i in xrange(n):
            running_sum += hidden_matrix[i][k] * X[i][j]
            denom += hidden_matrix[i][k]
        center.append(running_sum / denom)
    return center

# get centers based on hidden matrix
def get_centers(Data, hidden_matrix, k, n, m):
    X = []
    for j in xrange(k):
        X.append(get_center(Data, hidden_matrix, j, m))
    return X

# check if two sets of centers are equal
def are_centers_equal(X, old_X, k, m):
    for i in xrange(k):
        for j in xrange(m):
            if abs(X[i][j] - old_X[i][j]) > 0.0001:
                return False
    return True

# the fuzzy lloyd algorithm for kmeans
# Data: Data matrix
# k: number of centers
# beta: stiffness parameter
# n: number of rows in Data
# m: number of columns in Data
#
# returns: matrix of centers
def fuzzy_kmeans(Data, k, beta, n, m):
    maxtries = 1000

    old_X = [[0.] * m] * k    
    X = initialize(Data, k, n, m)
    i = 0
    
    while not are_centers_equal(X, old_X, k, m) and i < maxtries:
        i += 1
        old_X = X
        hidden_matrix = get_hidden_matrix(Data, X, beta, n, m)
        X = get_centers(Data, hidden_matrix, k, n, m)
    return X

# read problem data from file
# returns:
#   Data: Data matrix
#   k: number of centers
#   n: number of rows in Data
#   m: number of columns in Data
def readdat(filename):
    with open(filename, 'r') as f:
        Data = []
        k = int(f.readline().strip())
        beta = float(f.readline().strip())
        n = 0
        m = 0
    
        for string in f:
            curvec = string.strip().split()
            if m == 0:
                m = len(curvec)
            else:
                assert(m == len(curvec))
            Data.append(map(float, curvec))
            n += 1
    return((Data, k, beta, n, m))

# run program reading data from file
def main(filename):
    Data, k, beta, n, m = readdat(filename)
    centers = fuzzy_kmeans(Data, k, beta, n, m)
    print "\n".join([" ".join(map(lambda y: "%.1f" % y, x)) for x in centers])
    
if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
