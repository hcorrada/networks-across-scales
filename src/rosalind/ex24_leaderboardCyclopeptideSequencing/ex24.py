import sys
import itertools
import heapq

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

    def score(self, experimental_spectrum, cyclic=False):
        theoretical_spectrum = self.spectrum(cyclic=cyclic)

        score = 0
        
        for mass, count in theoretical_spectrum.iteritems():
            if mass in experimental_spectrum:
                if not cyclic:
                    if count <= experimental_spectrum[mass]:
                        score += count
                else:
                    score += min(count, experimental_spectrum[mass])
                     
        return score

def cut(leaderboard, n):
    if len(leaderboard) <= n:
        return
    
    cutoff = heapq.nlargest(n, leaderboard)[n-1]
    while leaderboard[0] < cutoff:
        heapq.heappop(leaderboard)

def expand(leaderboard, parent_mass, experimental_spectrum):
    out = []
    for score, peptide in leaderboard:
        for mass in masses:
            new_peptide = peptide.extend(mass)
            if new_peptide.total_mass() <= parent_mass:
                score = new_peptide.score(experimental_spectrum)
                heapq.heappush(out, (score, new_peptide))
    return out

def leaderboard_cyclopeptide_sequencing(spectrum, n):
    experimental_spectrum = {}
    for mass in spectrum:
        if mass in experimental_spectrum:
            experimental_spectrum[mass] += 1
        else:
            experimental_spectrum[mass] = 1
            
    parent_mass = max(spectrum)
    
    leaderboard = [(0, Peptide())]
    res = []
    best_peptide = None
    
    while len(leaderboard) > 0:
        leaderboard = expand(leaderboard, parent_mass, experimental_spectrum)
        for score, peptide in leaderboard:
            if peptide.total_mass() == parent_mass:
                cur_score = peptide.score(experimental_spectrum, cyclic=True)
                if best_peptide is None or cur_score > best_peptide[0]:
                    best_peptide = (cur_score, peptide)
        cut(leaderboard, n)
    return best_peptide[1]

filename = sys.argv[1]
with open(filename, 'r') as f:
    n = int(f.readline().strip())
    spectrum = map(int, f.readline().strip().split())
    peptide = leaderboard_cyclopeptide_sequencing(spectrum, n)
    print peptide
