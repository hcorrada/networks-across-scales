import sys
from Bio import SeqIO

# compute the z values for input string s
# z[i] = length of longest substring of s starting at position i
#        matching a prefix of s. I.E., s[:z[i]] = s[i:(i+z[i])]
# input:
#   s: a string
# output:
#   list of integer z values
def zalgorithm(s):
    n = len(s)
    z = [0] * n

    # start pointer of right-most z-box
    l = 0

    # find the end pointer of the right-most z-box
    # compare characters r to r+1 directly until the first
    # mismatch
    r = 0
    while r < n and s[r] == s[r+1]:
        r += 1

    # this gives z-value for second position
    z[1] = r

    # now compute the remaining z-values
    for k in xrange(2, n):
        # check for case 1, right-most z-box does not cover k
        if k > r:
            # compare characters directly until first mismatch
            t = 0
            while k+t < n and s[t] == s[k+t]:
                t += 1
            z[k] = t

            # this is the right-most z-box now,
            # so adjust start and end pointers accordingly
            if t > 0:
                l = k
                r = l + t - 1

        else:
            # case 2, right-most z-box does cover k, so k
            #   occurs within a substring that matches a prefix of s
            #
            # find the corresponding location to k in the prefix of s
            # and calculate the length of the remainder of matching prefix
            #   (this is string beta in Gusfield)
            kprime = k - l
            len_beta = r - k + 1

            # case 2a z[kprime] < |beta|
            # prefix match at kprime is shorter than beta
            # so prefix match at k must be the same length
            if z[kprime] < len_beta:
                z[k] = z[kprime]

            # case 2b z[kprime] > |beta|
            # prefix match at kprime is longer than beta
            # so prefix match at k is exactly beta
            elif z[kprime] > len_beta:
                z[k] = len_beta

            # case 2c z[kprime] == |beta|
            # now we are not sure about what the prefix match at k
            # looks like beyond beta, so we need to compare characters
            # and adjust the right-most z-box (which will now start at k)
            else:
                t = 1
                while r+t < n and s[r+t] == s[len_beta + t]:
                    t += 1
                z[k] = len_beta + t

                # adjust pointers for right-most z-box
                r += t - 1
                l = k
    return z

# compute the suffix-prefix-mismatch values from the
# z-values obtained via the fundamental pre-preprocessing algorithm
# spm[i]: length of longest substring of s ending at i that matches a prefix
#   of s (ending at j for example) AND s[i+1] != s[j+1]
# so, s[:spm[i]] = s[(i-spm[i]):i] AND s[spm[i]] != s[i+1]
def get_spm(z):
    n = len(z)

    # initialize spm
    spm = [0] * n

    # to fill-in spm[i] we need to find smallest value j
    # such that j + z[j] = i, so we'll do it by
    # going from the end of the string to the beginning
    # for each position j we find it's corresponding i and
    # set spm[i] = z[j], since we are moving down,
    # spm[i] may be overwritten by z[j] of a smaller j
    for j in reversed(xrange(2, n)):
        i = j + z[j] - 1
        spm[i] = z[j]

    return spm

# compute the KMP failure array from the spm values
# all we need to do is relax the requirement on mismatch character
# at the next position after spm[i], and we get the KMP failure function
def get_fa(spm):
    n = len(spm)
    fa = [0] * n

    # last value is given by spm
    fa[n-1] = spm[n-1]

    # now relax value from right to left
    for i in reversed(xrange(2, n-1)):
        fa[i] = max(fa[i+1] - 1, spm[i])
    return fa

def main(filename):
    s = SeqIO.read(filename, 'fasta').seq

    # use the fundamental pre-preprocessing algorithm
    # to compute z-values (this is the linear time preprocessing
    # algorithm from which we derive the KMP failure array)
    z = zalgorithm(s)

    # now convert z-values to suffix-prefix-mismatch values
    spm = get_spm(z)

    # finally convert spm values to failure array
    failure_array = get_fa(spm)
    print ' '.join(map(str, failure_array))


if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
