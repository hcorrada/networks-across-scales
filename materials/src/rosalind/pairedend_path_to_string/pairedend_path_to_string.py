import sys

# a class representing a node label
# using paired kmers
#
# slots:
#   _first: the first kmer in pair
#   _second: the second kmer in pair
class PairedKmer:
    # initialize from string with format
    # <first>|<second>
    def __init__(self, string):
        self._first, self._second = string.split("|")

    # return string representation of label with format
    # <first>|<second>
    def __repr__(self):
        return self._first + "|" + self._second

    # return the first kmer in pair
    def first(self):
        return self._first

    # return the second kmer in pair
    def second(self):
        return self._second

# a class to represent nodes in a DeBruijn graph
#
# slots:
#   _label: the paired k-1 mer label
#   _targets: list of target nodes for outgoing edges
class Node:
    # intialize node with label and empty targets list
    def __init__(self, label):
        self._label = label
        self._targets = []

    # get a string representation of node in format
    # <label> -> <comma_separated_list_of_target_labels>
    def __repr__(self):
        target_labels = [target.label_string() for target in self._targets]
        targets_string = ",".join(target_labels)
        return self.label_string() + " -> " + targets_string

    # return the node label
    def label(self):
        return self._label

    # return node label as string
    def label_string(self):
        return self.label().__repr__()

# class representing item in a doubly-linked list
# we use this for path representation as a doubly-linked
# list of nodes
class DoubleList:
    def __init__(self, node, prev_item=None, next_item=None):
        self._node = node
        self._prev_item = prev_item
        self._next_item = next_item

# class representing a path
#
# slots:
#   _k: length of each part of input paired k-mers
#   _d: gap in paired k-mers
#   _head: first item in doubly-linked list
#   _tail: last item in doubly_linked list
class Path:
    # initialize given k-mer size and gap
    # set head and tail to None
    def __init__(self, k, d):
        self._k = k
        self._d = d
        self._head = None
        self._tail = None

    # we use this helper method when reconstructing a string
    # from a given path. Given a list of kmers, constructs a string
    # by taking the first k-mer in the list and appending the last
    # character of remaining k-mers in succession
    #
    # input:
    #   kmers: a list of kmers
    #
    # output:
    #   string spelled out by kmer list
    @staticmethod
    def _string_helper(kmers):
        out = ""
        for kmer in kmers:
            # use the complete kmer if it's the first
            # i.e., when output string is empty
            if len(out) == 0:
                out = kmer
            else:
                # otherwise append last character of kmer
                out += kmer[-1]
        return out

    # get string spelled out by path of paired k-mers
    # output:
    #   string spelled out by paired-kmers in path
    def get_string(self):
        # get string spelled out by first k-mer in each pair
        prefix_string = self._string_helper([node.label().first() for node in self.nodes()])

        # get string spelled out by second k-mer in each pair
        suffix_string = self._string_helper([node.label().second() for node in self.nodes()])

        # check the overlap between prefix and suffix strings
        # is valid
        n = len(prefix_string)
        prefix_overlap = prefix_string[self._k + self._d:]
        suffix_overlap = suffix_string[:n - self._k - self._d]

        if prefix_overlap != suffix_overlap:
            # invalid overlap, no string is spelled out by this path
            return None
        else:
            # valid overlap, return the string
            return prefix_string + suffix_string[-(self._k + self._d):]

    # return true if path is empty
    def is_empty(self):
        return self._head is None and self._tail is None

    # append node to the tail of the list
    def append(self, node):
        # create a new item
        item = DoubleList(node)
        if self.is_empty():
            # make new item head and tail of list
            self._head = self._tail = item
        else:
            # grab the item in the tail
            tail_item = self._tail

            # make item in tail point to new item
            tail_item._next_item = item

            # make new item's back pointer point to item in tail
            item._prev_item = tail_item

            # make the new item the tail of the list
            self._tail = item

    # a generator for nodes in the path
    def nodes(self):
        current = self._head
        while current is not None:
            yield current._node
            current = current._next_item

    # a string representation of the path
    def __repr__(self):
        node_labels = [node._label for node in self.nodes()]
        return " -> ".join([label.__repr__() for label in node_labels])

# make a Path object from a list of edges
def build_pairedend_path(paired_kmers, k, d):
    # initialize path object
    path = Path(k, d)

    # create paired kmer objects from string
    paired_kmers = map(PairedKmer, paired_kmers)

    # append paired kmer objects to path object
    for paired_kmer in paired_kmers:
        node = Node(paired_kmer)
        path.append(node)

    return path

# read input from file
def readdat(filename):
    with open(filename, 'r') as f:
        k, d = [int(s) for s in f.readline().strip().split()]
        paired_kmers = []
        for line in f:
            paired_kmers.append(line.strip())
        return paired_kmers, k, d

def main(filename):
    paired_kmers, k, d = readdat(filename)
    path = build_pairedend_path(paired_kmers, k, d)
    string = path.get_string()
    print string

# this is here so this plays nicely with ipython %loadpy magic
if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
