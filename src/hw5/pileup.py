from collections import defaultdict
import math

class PileUpStrategy:
    IGNORE = 0
    UNIFORM = 1
    FIRST = 2
    

class PileUp:
    def __init__(self, text, strategy=PileUpStrategy.IGNORE):
        self._reference = text
        self._strategy = strategy
        self._pileup = [defaultdict(int) for i in xrange(len(text))]
        
    def insert(self, positions, read):
        if len(positions) > 1 and self._strategy == PileUpStrategy.IGNORE:
            return

        position = positions[0]

        
        if len(positions) > 1 and self._strategy == PileUpStrategy.UNIFORM:
            position = random.choice(positions)

        for c in list(read):
            self._pileup[position][c] += 1
            position += 1

    def __str__(self):
        out = ''
        
        i = -1
        for d in self._pileup:
            i += 1
            if len(d) == 0:
                continue
            if any([nuc != self._reference[i] for nuc in d.iterkeys()]):
                out += '*'
                
            out += "%d %s " % (i, self._reference[i])
            calls = ["%s:%d" % (nuc,count) for nuc, count in d.iteritems()]
            out +=  ','.join(calls) + '\n'
        return out

    def mismatching_positions(self, minrate=0.0):
        out = []
        i = -1
        for d in self._pileup:
            i += 1
            if len(d) == 0:
                continue
            total = sum(d.values())
            cutoff = math.floor(total * minrate)
            mismatches = [(nuc, count) for nuc,count in d.iteritems() if nuc != self._reference[i] and count >= cutoff]
            if len(mismatches) > 0:
                out.append((i, self._reference[i], mismatches))
        return out

    def print_mismatches(self, minrate=0.0):
        positions = self.mismatching_positions(minrate)
        for p, nuc, mismatches in positions:
            print p, nuc,
            print ','.join(["%s:%d" % mismatch for mismatch in mismatches])
        return
            
    
            
