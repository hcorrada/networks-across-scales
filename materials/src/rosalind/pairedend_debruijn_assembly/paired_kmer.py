# a class representing a node label
# using paired kmers
#
# slots:
#   _first: the first kmer in pair
#   _second: the second kmer in pair
class PairedKmer:
    # initialize from string with format
    # <first>|<second>
    def __init__(self, string, second=None):
        if second is None:
            self._first, self._second = string.split("|")
        else:
            self._first = string
            self._second = second

    def as_string(self):
        return self._first + "|" + self._second

    # return string representation of label with format
    # <first>|<second>
    def __repr__(self):
        return self.as_string()

    # how to hash this, we need this to
    # use this as dictionary key
    def __hash__(self):
        return hash(self.as_string())

    # equality based on string equality
    def __eq__(self, other):
        return self.as_string() == other.as_string()

    # return the first kmer in pair
    def first(self):
        return self._first

    # return the second kmer in pair
    def second(self):
        return self._second

    # return the paired k-1 mer prefix
    def prefix(self):
        return PairedKmer(self._first[:-1], self._second[:-1])

    # return the paired k-1 mer suffix
    def suffix(self):
        return PairedKmer(self._first[1:], self._second[1:])
