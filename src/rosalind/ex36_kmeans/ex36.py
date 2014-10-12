import sys
import heapq
import math

# compute euclidean distance between
# vectors x and y of length m
def point_to_point_distance(x, y, m):
    res = 0
    for i in xrange(m):
        res += math.pow(x[i] - y[i], 2)
    return res

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

# assign point to it's closest center in X
def assign_point(y, X, n, m):
    heap = []
    for i in xrange(len(X)):
        heapq.heappush(heap, (point_to_point_distance(y, X[i], m), i))

    # smallest distance in root of heap with entry (distance, center index)
    return heap[0][1]

# assign each point in Data to it's closest center in X
def assign(Data, X, n, m):
    clusters = []
    for i in xrange(n):
        clusters.append(assign_point(Data[i], X, n, m))
    return clusters

# get center of a set of points
def get_center(X, m):
    center = []
    n = len(X)
    
    for j in xrange(m):
        sum = 0.
        for i in xrange(n):
            sum += X[i][j]
        center.append(sum / n)
    return center

# get centers based on clusters
def get_centers(Data, clusters, k, n, m):
    X = []
    for i in xrange(k):
        Y = []
        for j in xrange(n):
            if clusters[j] == i:
                Y.append(Data[j])
        X.append(get_center(Y, m))
    return X

# check if two cluster assignmetns are equal
def are_clusters_equal(x, y):
    for i in xrange(len(x)):
        if x[i] != y[i]:
            return False
    return True

# the lloyd algorithm for kmeans
# Data: Data matrix
# k: number of centers
# n: number of rows in Data
# m: number of columns in Data
#
# returns: matrix of centers
def kmeans(Data, k, n, m):
    maxtries = 1000
    
    X = initialize(Data, k, n, m)
    clusters = assign(Data, X, n, m)
    old_clusters = [-1] * n
    i = 0
    
    while not are_clusters_equal(clusters, old_clusters) and i < maxtries:
        i += 1
        old_clusters = clusters
        X = get_centers(Data, clusters, k, n, m)
        clusters = assign(Data, X, n, m)
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
    return((Data, k, n, m))

# run program reading data from file
def main(filename):
    Data, k, n, m = readdat(filename)
    centers = kmeans(Data, k, n, m)
    print "\n".join([" ".join(map(lambda y: "%.2f" % y, x)) for x in centers])
    
if __name__ == '__main__':
    filename = sys.argv[1]
    main(filename)
