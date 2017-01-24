import numpy
import sys

def printout(S, k, l, r, z):
    print 'k=%d, l=%d, r=%d' % (k,l,r)
    print S
    print ''.join(map(lambda x: str(int(x)), z))

    sys.stdout.write(''.join([' ' for i in xrange(l)]))
    sys.stdout.write('l')

    if r>l+1:
        sys.stdout.write(''.join([' ' for i in xrange(l+1,r)]))
    print 'r'

    sys.stdout.write(''.join([' ' for i in xrange(k)]))
    print 'k'

    print '\n'
    
def zScores(S):
    n = len(S)
    z = numpy.zeros([n])

    l = 0


    # compute z2
    k = 1
    r = k
    printout(S, k, l, r, z)

    while S[r] == S[r-1] and r < n:
        print '%s=%s' % (S[r], S[r-1])
        z[1] += 1
        r += 1

    # decrease r so it points to last matching position
    if r > 1:
        r -= 1    
        l = k

        
    for k in xrange(2,n):
        printout(S, k, l, r, z)
        
        if r < k:
            # case 1 there is no z-box covering k
            t = 0
            while S[t] == S[k+t] and k+t < n:
                t += 1

            z[k] = t

            if t > 0:
                r = k + t - 1
                l = k

        else:
            kprime = k - l
            beta = r - k + 1
            if z[kprime] < beta:
                z[k] = z[kprime]
            elif z[kprime] > beta:
                z[k] = beta
            else:
                t = 1
                while S[r+t] == S[beta + t] and r+t < n:
                    t += 1
                z[k] = beta + t -1
                if t > 1:
                    r += t - 1
                    l = k
    printout(S,k,l,r,z) 
    return z

def main(S):
    print zScores(S)

if __name__ == '__main__':
    main('ABCDABD$ABCABCDABABCDABDABDE')


