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
    def __init__(self):
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
    def get_string(self, k, d):
        # get string spelled out by first k-mer in each pair
        prefix_string = self._string_helper([node.label().first() for node in self.nodes()])

        # get string spelled out by second k-mer in each pair
        suffix_string = self._string_helper([node.label().second() for node in self.nodes()])

        # check the overlap between prefix and suffix strings
        # is valid
        n = len(prefix_string)
        prefix_overlap = prefix_string[k + d:]
        suffix_overlap = suffix_string[:-(k + d)]

        if prefix_overlap != suffix_overlap:
            # invalid overlap, no string is spelled out by this path
            return None
        else:
            # valid overlap, return the string
            return prefix_string + suffix_string[-(k + d):]

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

    def find_item(self, node):
        cur_item = self._head
        while cur_item is not None:
            if cur_item._node.label() == node.label():
                return cur_item
            cur_item = cur_item._next_item
        return None

    def stitch(self, item, other):
        if self.is_empty():
            return other
        if item is None:
            return self

        # get reference to previous item on this path
        prev_item = item._prev_item

        # stick head of other path after previous item on this path
        prev_item._next_item = other._head
        other._head._prev_item = prev_item

        # put tail of other path before this item
        other._tail._next_item = item
        item._prev_item = other._tail

        return self

    # a generator for nodes in the path
    def nodes(self):
        current = self._head
        while current is not None:
            yield current._node
            current = current._next_item
            if current == self._head:
                current = None

    # a string representation of the path
    def __repr__(self):
        node_labels = [node._label for node in self.nodes()]
        return " -> ".join([label.__repr__() for label in node_labels])
