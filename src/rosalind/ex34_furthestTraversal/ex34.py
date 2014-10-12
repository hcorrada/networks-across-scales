import sys
import math
import heapq

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

# the furthest traversal algorithm
# Data: data matrix
# k: number of centers to return
# n: number of rows in Data
# m: number of columns in Data
def furthest_traversal(Data, k, n, m):
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


# run the program
def main(filename):
    Data, k, n, m = readdat(filename)
    centers = furthest_traversal(Data, k, n, m)
    print "\n".join([" ".join(map(str,x)) for x in centers])

if __name__ == '__main__' and "get_ipython" not in dir():
    filename=sys.argv[1]
    main(filename)

    
    
                
