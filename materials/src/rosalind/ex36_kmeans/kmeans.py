import sys
import heapq
import math

# compute squared euclidean distance between
# vectors x and y of length m
def point_to_point_distance(x, y, m):
    res = 0.
    for i in xrange(m):
        res += math.pow(x[i] - y[i], 2)
    return res

# initialize the set of centers
# by selecting the first k points in data
def initialize_centers(data, k, n, m):
    centers = []
    for i in xrange(k):
        centers.append(data[i])
    return centers

# assign point to its closest center in centers
def assign_point(y, centers, k, m):
    heap = []
    for i in xrange(k):
        heapq.heappush(heap, (point_to_point_distance(y, centers[i], m), i))

    # smallest distance in root of heap with entry (distance, center index)
    return heap[0][1]

# assign each point in data to its closest center in centers
def centers_to_clusters(data, centers, k, n, m):
    clusters = []
    for i in xrange(n):
        clusters.append(assign_point(data[i], centers, k, m))
    return clusters

# compute center of a set of points
def compute_center(points, m):
    center = []
    n = len(points)

    for j in xrange(m):
        cur_sum = 0.
        for i in xrange(n):
            cur_sum += points[i][j]
        center.append(cur_sum / n)
    return center

# get centers based on clusters
def clusters_to_centers(data, clusters, k, n, m):
    centers = []
    for i in xrange(k):
        points_in_cluster = []
        for j in xrange(n):
            if clusters[j] == i:
                points_in_cluster.append(data[j])
        new_center = compute_center(points_in_cluster, m)
        centers.append(new_center)
    return centers

# check if two cluster assignments are equal
def are_clusters_equal(x, y, n):
    for i in xrange(n):
        if x[i] != y[i]:
            return False
    return True

# the lloyd algorithm for kmeans
# data: Data matrix
# k: number of centers
# n: number of rows in Data
# m: number of columns in Data
#
# returns: matrix of centers
def kmeans(data, k, n, m):
    maxtries = 1000

    centers = initialize_centers(data, k, n, m)
    clusters = centers_to_clusters(data, centers, k, n, m)
    old_clusters = [-1] * n
    i = 0

    while not are_clusters_equal(clusters, old_clusters, n) and i < maxtries:
        i += 1
        old_clusters = clusters
        centers = clusters_to_centers(data, clusters, k, n, m)
        clusters = centers_to_clusters(data, centers, k, n, m)
    return centers

# read problem data from file
# returns:
#   data: Data matrix
#   k: number of centers
#   n: number of rows in Data
#   m: number of columns in Data
def readdat(filename):
    with open(filename, 'r') as f:
        k, m = map(int, f.readline().strip().split())
        n = 0
        data = []

        for line in f:
            curvec = map(float, line.strip().split())
            assert(len(curvec) == m)
            data.append(curvec)
            n += 1
    return data, k, n, m

# run program reading data from file
def main(filename):
    data, k, n, m = readdat(filename)
    centers = kmeans(data, k, n, m)
    center_to_string = lambda center: " ".join(map(lambda x: "%.3f" % x, center))
    print "\n".join(map(center_to_string, centers))

if __name__ == '__main__':
    filename = sys.argv[1]
    main(filename)
