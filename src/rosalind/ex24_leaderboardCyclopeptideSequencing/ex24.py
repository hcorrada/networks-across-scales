import sys
import itertools
import heapq

# the list of (unique) amino acid masses 
masses = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]

# class representing a peptide
#   has a single slot '_chain' containing the list of masses corresponding to the peptide
class Peptide:
    # constructor: the empty peptide if mass is None, otherwise a single aminoacid peptide
    def __init__(self, mass=None):
        if mass is None:
            self._chain = []
        else:
            self._chain = [mass]

    # make a nice string representing the peptide
    def __repr__(self):
        return "-".join(map(str, self._chain))

    # make a new peptide by extending this peptide with given mass
    def extend(self, mass):
        # make a new peptide
        res = Peptide()
        # store a copy of the aminoacid chain for this peptide
        # in the new peptide
        res._chain = list(self._chain)
        # appen the new mass to the chain of the new peptide
        res._chain.append(mass)
        return res

    # the length of the this peptide
    def __len__(self):
        return len(self._chain)

    # the total mass of this peptide: the sum of masses
    def total_mass(self):
        return sum(self._chain)

    # return the spectrum for this peptide
    # if cyclic == True, then return the cyclospectrum, otherwise return the linear spectrum
    # we store spectrums as hash tables: keys are masses, values are number of times that mass occurs
    # in the spectrum
    def spectrum(self, cyclic=False):
        k = len(self)
        # initialize spectrum with 0 mass
        res = {0: 1}

        # get all possible start,end points of subpeptides for this peptide
        for i,j in itertools.combinations(xrange(k+1), 2):
            # compute mass for subpeptide
            mass = sum(self._chain[i:j])
            # add to hash table if needed, otherwise increase count
            if mass in res:
                res[mass] += 1
            else:
                res[mass] = 1

        # done if not cyclic
        if not cyclic:
            return res

        # if cyclic, add masses for cyclic version of the peptide
        for i in xrange(1, k):
            a = sum(self._chain[i:])
            for j in xrange(i-1):
                mass = a + sum(self._chain[:j+1])
                if mass in res:
                    res[mass] += 1
                else:
                    res[mass] = 1
        return res

    # score this peptide against the given experimental spectrum
    # input:
    #   experimental_spectrum: hash table
    #   cyclic: boolean, compute as cyclic?
    def score(self, experimental_spectrum, cyclic=False):
        # get the theoretical spectrum for this peptide
        theoretical_spectrum = self.spectrum(cyclic=cyclic)

        score = 0

        # go through masses in the spectrum for this peptide        
        for mass, count in theoretical_spectrum.iteritems():
            # check if this mass is in the experimental spectrum
            if mass in experimental_spectrum:
                # it is, we score depending on cyclic or not
                if not cyclic:
                    # if it is not cyclic, increase score as long as the count for this peptide
                    # is less than or equal to the count in the experimental spectrum
                    if count <= experimental_spectrum[mass]:
                        score += count
                else:
                    # if scoring as cyclic, score equals the smaller of observed and theoretical counts
                    score += min(count, experimental_spectrum[mass])
                     
        return score

# keep only peptides with top n scores in the leaderboard
# input:
#   leaderboard: a heap
#   n: the top score to keep
def cut(leaderboard, n):
    # don't do anything is leaderboard is too small
    if len(leaderboard) <= n:
        return

    # find out the top n-th score in leaderboard
    cutoff = heapq.nlargest(n, leaderboard)[n-1]
    # pop root of heap as long as it's smaller than the cutoff
    while leaderboard[0] < cutoff:
        heapq.heappop(leaderboard)

# expand each peptide in the leaderboard
# input:
#  leaderboard: a heap
#  parent_mass: the maximum allowed mass for a peptide
#  experimental_spectrum: the observered spectrum
def expand(leaderboard, parent_mass, experimental_spectrum):
    # a new heap
    out = []
    # expand each peptide in the leaderboard
    for score, peptide in leaderboard:
        for mass in masses:
            # make the new peptide
            new_peptide = peptide.extend(mass)
            # if it's not too heavy
            if new_peptide.total_mass() <= parent_mass:
                # score the (linear) peptide 
                score = new_peptide.score(experimental_spectrum)
                # add to the new heap
                heapq.heappush(out, (score, new_peptide))
    return out

# heuristic, greedy cyclopeptide sequencing algorithm
# keeps best scoring n peptides, at each peptide length
# input:
#  spectrum: observed experimental spectrum
#  n: size of leaderboard
# output:
#  best scoring peptide
def leaderboard_cyclopeptide_sequencing(spectrum, n):
    # store the experimental spectrum as a hash table
    # keys are masses, values are number of times each mass occurs in spectrum
    experimental_spectrum = {}
    for mass in spectrum:
        if mass in experimental_spectrum:
            experimental_spectrum[mass] += 1
        else:
            experimental_spectrum[mass] = 1

    # get the mass of full peptide
    parent_mass = max(spectrum)

    # initialize heap (leaderboard) with empty peptide
    leaderboard = [(0, Peptide())]
    res = []
    best_peptide = None

    # while the leaderboard is not empty
    while len(leaderboard) > 0:
        # make a new leaderboard by expanding (adding one mass) to the
        # peptides in the current leaderboard
        # only adds new peptides if they are not too heavy
        # adds new peptides into heap based on score vs. experimental spectrum
        leaderboard = expand(leaderboard, parent_mass, experimental_spectrum)

        # check if the current leaderboard contains a candidate solution
        for score, peptide in leaderboard:
            # ok, this peptide matches the parent mass, so it's a candidate solution
            if peptide.total_mass() == parent_mass:
                # score the peptide, now we do it as a cyclic peptide since we know it has the right parent mass
                cur_score = peptide.score(experimental_spectrum, cyclic=True)
                # check if this score is better than current best
                if best_peptide is None or cur_score > best_peptide[0]:
                    best_peptide = (cur_score, peptide)
        # cut the leaderboard to top n scoring peptides
        cut(leaderboard, n)
    return best_peptide[1]

filename = sys.argv[1]
with open(filename, 'r') as f:
    n = int(f.readline().strip())
    spectrum = map(int, f.readline().strip().split())
    peptide = leaderboard_cyclopeptide_sequencing(spectrum, n)
    print peptide
