import sys
import itertools
from collections import deque

masses = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]

class Peptide:
    def __init__(self, mass=None):
        if mass is None:
            self._chain = []
        else:
            self._chain = [mass]

    def __repr__(self):
        return "-".join(map(str, self._chain))

    def extend(self, mass):
        res = Peptide()
        res._chain = list(self._chain)
        res._chain.append(mass)
        return res

    def __len__(self):
        return len(self._chain)

    def total_mass(self):
        return sum(self._chain)
    
    def spectrum(self, cyclic=False):
        k = len(self)
        res = {0: 1}

        for i,j in itertools.combinations(xrange(k+1), 2):
            mass = sum(self._chain[i:j])
            if mass in res:
                res[mass] += 1
            else:
                res[mass] = 1

        if not cyclic:
            return res

        for i in xrange(1, k):
            a = sum(self._chain[i:])
            for j in xrange(i-1):
                mass = a + sum(self._chain[:j+1])
                if mass in res:
                    res[mass] += 1
                else:
                    res[mass] = 1
        return res

    def consistent(self, experimental_spectrum):
        theoretical_spectrum = self.spectrum()
        
        for mass in theoretical_spectrum.keys():
            if mass not in experimental_spectrum:
                return False
            if theoretical_spectrum[mass] > experimental_spectrum[mass]:
                return False
                
        return True
    
def cyclopeptide_sequencing(spectrum):
    experimental_spectrum = {}
    for mass in spectrum:
        if mass in experimental_spectrum:
            experimental_spectrum[mass] += 1
        else:
            experimental_spectrum[mass] = 1
            
    parent_mass = max(spectrum)
    
    q = deque()
    for mass in masses:
        a = Peptide(mass)
        if a.consistent(experimental_spectrum):
            q.append(a)

    res = []
    while len(q) > 0:
        cur_peptide = q.popleft()
        for mass in masses:
            new_peptide = cur_peptide.extend(mass)
            
            if new_peptide.total_mass() == parent_mass and new_peptide.spectrum(cyclic=True) == experimental_spectrum:
                 res.append(new_peptide)
            elif new_peptide.consistent(experimental_spectrum):
                q.append(new_peptide)
            
    return res

filename = sys.argv[1]
with open(filename, 'r') as f:
    spectrum = map(int, f.read().strip().split())
    peptides = cyclopeptide_sequencing(spectrum)
    print " ".join(map(str, peptides))
