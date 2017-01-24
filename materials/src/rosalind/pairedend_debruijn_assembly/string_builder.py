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
def build_string(path, k, d):
    # get string spelled out by first k-mer in each pair
    first_kmers = [node.label().first() for node in path.nodes()]
    prefix_string = _string_helper(first_kmers)

    # get string spelled out by second k-mer in each pair
    second_kmers = [node.label().second() for node in path.nodes()]
    suffix_string = _string_helper(second_kmers)

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
