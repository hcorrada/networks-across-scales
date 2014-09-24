import sys
import itertools

keys = list("GASPVTCILNDKQEMHFRYW")
masses = [57, 71, 87, 97, 99, 101, 103, 113, 113, 114, 115, 128, 128, 129, 131, 137, 147, 156, 163, 186]
mass_table = dict(zip(keys, masses))

def weigh_it(peptide):
    return sum([mass_table[x] for x in peptide])

def theoretical_spectrum(peptide):
    k = len(peptide)
    res = [0]

    # compute weights for all peptides in the linear version    
    for i,j in itertools.combinations(xrange(k+1), 2):
        res.append(weigh_it(peptide[i:j]))

    # compute for circularized versions
    # for each starting position
    for i in xrange(1, k):
        # take suffix starting at that position
        a = peptide[i:]

        # take all valid substrings from 1:that position - 1
        for j in xrange(i-1):
            b = peptide[:j+1]
            res.append(weigh_it(a+b))
    return sorted(res)

filename = sys.argv[1]
with open(filename, 'r') as f:
    peptide = f.readline().strip()
    print " ".join(map(str, theoretical_spectrum(peptide)))
