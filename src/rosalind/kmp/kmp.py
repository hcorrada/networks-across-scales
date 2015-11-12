import sys
from Bio import SeqIO

def zalgorithm(s):
    n = len(s)
    z = [0] * n

    l = 0

    # compute z2
    k = 1
    r = k-1

    while s[r] == s[r+1] and r < n:
        z[1] += 1
        r += 1

    for k in xrange(2, n):
        # check for case 1, no z box covering k
        if k > r:
            # compare strings directly
            t = 0
            while s[t] == s[k+t] and k+t < n:
                t += 1
            z[k] = t

            if t > 0:
                l = k
                r = l + t - 1

        else:
            # case 2, rightmost zbox covers k
            kprime = k - l
            len_beta = r - k + 1

            # case 2a z[kprime] < |beta|
            if z[kprime] < len_beta:
                z[k] = z[kprime]

            # case 2b z[kprime] > |beta|
            elif z[kprime] > len_beta:
                z[k] = len_beta

            # case 2c z[kprime] == |beta|
            else:
                t = 1
                while s[r+t] == s[len_beta + t] and r+t < n:
                    t += 1
                z[k] = len_beta + t
                r += t - 1
                l = k
    return z

def get_failure_func(z):
    n = len(z)
    sp_prime = [0] * n

    for k in xrange(2, n):
        j = n - k + 1
        i = j + z[j] - 1
        sp_prime[i] = z[j]

    spm = sp_prime
    for k in xrange(3, n):
        i = n - k + 1
        spm[i] = max(spm[i+1] - 1, sp_prime[i])
    return spm

def main(filename):
    s = SeqIO.read(filename, 'fasta').seq
    z = zalgorithm(s)
    failure_func = get_failure_func(z)
    print ' '.join(map(str, failure_func))


if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
