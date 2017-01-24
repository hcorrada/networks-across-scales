import sys
from collections import deque

masses = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]

def count_peptides(m):
    active_set = {}
    total_count = 0
        
    for mass in masses:
        if mass < m:
            active_set[mass] = 1
        if mass == m:
            total_count += 1

    while len(active_set) > 0:
        cur_mass = min(active_set)
        count = active_set.pop(cur_mass)
        
        for mass in masses:
            s = cur_mass + mass
            
            if s == m:
                total_count += count
                
            if s < m:
                if s in active_set:
                    active_set[s] += count
                else:
                    active_set[s] = count
    return total_count

filename = sys.argv[1]
with open(filename, 'r') as f:
    m = int(f.readline().strip())
    print count_peptides(m)
