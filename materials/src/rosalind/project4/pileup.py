from collections import defaultdict
import math

# flag to indicate how to handle reads that map to multiple positions
# in the reference
class PileUpStrategy:
    IGNORE = 0 # ignore these reads
    UNIFORM = 1 # select a position uniformly at random
    FIRST = 2 # sleect the first position in list


# class to print pileups of aligned reads to a given reference
#
# Fields:
#   _reference: the reference string
#   _strategy: how to deal with reads that map to multiple positions
#       see class 'PileUpStrategy'
#   _pileup: the pileup itself, an array of dictionaries, one per position
#       in the reference, each dictionary holds the number of times each
#       nucleotide has been aligned to that position in a set of reads
class PileUp:
    def __init__(self, text, strategy=PileUpStrategy.IGNORE):
        self._reference = text
        self._strategy = strategy

        # initialize array of dictionaries
        self._pileup = [defaultdict(int) for i in xrange(len(text))]

    # update the pileup given a set of alignments
    #
    # Input:
    #   positions: indices to which read has been aligned
    #   read: the read sequence
    #
    # reads that map to multiple positions are handled according
    # to self._strategy
    def insert(self, positions, read):
        # more than one position, and we are supposed to ignore
        # reads that map to more than one position
        if len(positions) > 1 and self._strategy == PileUpStrategy.IGNORE:
            return

        # take the first position to which the read mapped
        position = positions[0]

        # we are supposed to choose between multiple postiions at random
        # and we have more than one position
        if self._strategy == PileUpStrategy.UNIFORM and len(positions) > 1:
            position = random.choice(positions)

        # for each character in the read, update the
        # dictionary at the respective position
        for c in read:
            self._pileup[position][c] += 1
            position += 1

    # print pileup representation
    def __str__(self):
        out = ''

        i = -1
        # for each position in reference
        for d in self._pileup:
            i += 1

            # no bases have aligned here
            if len(d) == 0:
                continue

            # there is at least one mismatch here
            if any([nuc != self._reference[i] for nuc in d.iterkeys()]):
                out += '*'

            # print out position, reference character
            out += "%d %s " % (i, self._reference[i])

            # and then nucleotides observed in reads and the number of times
            # observed
            calls = ["%s:%d" % (nuc,count) for nuc, count in d.iteritems()]
            out +=  ','.join(calls) + '\n'
        return out

    # find positions in pileup that have mismatches occurring at
    # frequency greater than given rate
    #
    # Input:
    #   minrate: float, minimum mismatch frequency required to include position
    # Output:
    #   array of tuples [(pos, ref_nuc, [(read_nuc, count)])])]:
    #       pos: position in reference
    #       ref_nuc: reference nucleotide at given position
    #       read_nuc: mismatching nucleotide in read
    #       count: number of times mismatching nucleotide is observed
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

    # print mismatching positions in pileup occuring at
    # frequency greater than given rate
    #
    # Input:
    #   minrate: float, minimum mismatch frequency required to print
    def print_mismatches(self, minrate=0.0):
        # find positions with mismatches
        positions = self.mismatching_positions(minrate)
        # now print them out
        for p, nuc, mismatches in positions:
            print p, nuc,
            print ','.join(["%s:%d" % mismatch for mismatch in mismatches])
        return
