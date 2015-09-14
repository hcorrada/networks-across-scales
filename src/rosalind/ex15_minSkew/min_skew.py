import sys

# compute minimum skew for given string
# input:
#  genome: string
#
# output:
#  the set of indices with minimum skew
#
# this is a linear algorithm: adjusts running skew in each position
def minimum_skew(genome):
    running_skew = 0
    cur_min = len(genome)+1
    best_indices = [0]

    for i in xrange(len(genome)):
        # take the current character
        c = genome[i]

        # decrease skew if C
        if c == 'C':
            running_skew = running_skew - 1
        # increase skew if G
        elif c == 'G':
            running_skew = running_skew + 1

        # add this position if equal to current minimum
        if running_skew == cur_min:
            best_indices.append(i + 1)
        # this has smaller skew, start new list with this position
        elif running_skew < cur_min:
            cur_min = running_skew
            best_indices = [i + 1]

    return best_indices
            
filename = sys.argv[1]
with open(filename, 'r') as f:
    genome = f.read().strip()
    print " ".join(map(str, minimum_skew(genome)))
