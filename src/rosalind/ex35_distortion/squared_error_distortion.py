import sys
import re
import heapq
import math

# compute squared euclidean distance between
# vectors x and y of length m
def point_to_point_distance(x, y, m):
    res = 0.
    for i in xrange(m):
        res += math.pow(x[i] - y[i], 2)
    return res

# return the minimum euclidean distance
# between point y and set of k points centers
# points are of dimension m
def point_to_set_distance(y, centers, k, m):
    heap = []
    for i in xrange(k):
        heapq.heappush(heap, point_to_point_distance(y, centers[i], m))
    return heap[0]

# compute squared error distortion
# average distance between points in data to closest center in centers
def squared_error_distortion(data, centers, nrows, k, m):
    res = 0.
    for i in xrange(nrows):
        res += point_to_set_distance(data[i], centers, k, m)

    return res / nrows

# read problem data from file
# returns:
#   data: Data matrix
#   centers: center matrix
#   nrows: number of rows in Data
#   k: number of centers
#   m: number of data dimensions
def readdat(filename):
    with open(filename, 'r') as f:
        k, m = map(int, f.readline().strip().split())
        centers = []
        for _ in xrange(k):
            cur_center = map(float, f.readline().strip().split())
            centers.append(cur_center)

        # skip this line
        _ = f.readline()

        data = []
        nrows = 0

        # now read the remaining lines
        for line in f:
            nrows += 1
            cur_point = map(float, line.strip().split())
            data.append(cur_point)

        return data, centers, nrows, k, m

def main(filename):
    data, centers, nrows, k, m = readdat(filename)
    dist = squared_error_distortion(data, centers, nrows, k, m)
    print "%.3f" % dist

if __name__ == '__main__':
    filename = sys.argv[1]
    main(filename)
