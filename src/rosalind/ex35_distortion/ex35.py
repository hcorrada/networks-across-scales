import sys
import re
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

# compute squared error distortion
# average distance between points in Data to closest center in X
def squared_error_distortion(Data, X, ndata, nx, m):
    res = 0.
    for i in xrange(ndata):
        res += point_to_set_distance(Data[i], X, m)
        
    return res / ndata

# read problem data from file
# returns:
#   Data: Data matrix 
#   X: center matrix
#   ndata: number of rows in Data
#   nx: number of rows in X
#   m: number of columns in Data and X
def readdat(filename):
    X = []
    Data = []

    ndata = 0
    nx = 0
    m = 0

    reading_x = True
    
    with open(filename, 'r') as f:
        for string in f:
            string = string.strip()
            if re.match("^Dat", string) is not None:
                reading_x = False
                continue

            vec = string.split()
            if m == 0:
                m = len(vec)
            else:
                assert(len(vec) == m)

            if reading_x:
                X.append(map(float, vec))
                nx += 1
            else:
                Data.append(map(float, vec))
                ndata += 1        
    return (Data, X, ndata, nx, m)

def main(filename):
    Data, X, ndata, nx, m = readdat(filename)
    dist = squared_error_distortion(Data, X, ndata, nx, m)
    print "%.3f" % dist
    
if __name__ == '__main__':
    filename = sys.argv[1]
    main(filename)
