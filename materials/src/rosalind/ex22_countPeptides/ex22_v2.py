import sys
from collections import deque

masses = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]

tab = {}

def helper(m):
    if m < 57:
        return 0
    
    count = 1 if m in masses else 0
    for mass in masses:
        if m-mass in tab:
            count += tab[m - mass]
        else:
            count += helper(m-mass)
            
    tab[m] = count
    return tab[m]

def count_peptides(m):
    return helper(m)
    
filename = sys.argv[1]
with open(filename, 'r') as f:
    m = int(f.readline().strip())
    print count_peptides(m)
