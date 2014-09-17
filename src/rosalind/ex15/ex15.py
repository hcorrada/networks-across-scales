import sys

def minimum_skew(genome):
    running_skew = 0
    cur_min = len(genome)+1
    best_indices = [0]

    for i in xrange(len(genome)):
        c = genome[i]
        
        if c == 'C':
            running_skew = running_skew - 1
        elif c == 'G':
            running_skew = running_skew + 1

        if running_skew == cur_min:
            best_indices.append(i + 1)
        elif running_skew < cur_min:
            cur_min = running_skew
            best_indices = [i + 1]

    return best_indices
            
filename = sys.argv[1]
with open(filename, 'r') as f:
    genome = f.read().strip()
    print " ".join(map(str, minimum_skew(genome)))
