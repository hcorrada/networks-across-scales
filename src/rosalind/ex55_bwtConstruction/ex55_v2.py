import sys

def make_cmp(text):
    n = len(text)
    
    def _cmp(x,y):
        d = cmp(text[x], text[y])
        i = x
        j = y
        k = 0
        while d == 0 and k < n:
            i = (i + 1) % n
            j = (j + 1) % n

            d = cmp(text[i], text[j])
        return d
    return _cmp
                
# construct the burrows-wheeler transform of text
def get_bwt(text):
    # setup rotation index array
    indices = range(len(text))

    # sort index array by string rotations
    # uses comparison function defined above
    indices = sorted(indices, cmp=make_cmp(text))

    # get the last character of each rotation
    bwt = [text[i-1] for i in indices]

    # turn into a string and return
    return ''.join(bwt)

def readdat(filename):
    with open(filename, 'r') as f:
        return f.readline().strip()

def main(filename):
    text = readdat(filename)
    bwt = get_bwt(text)
    print bwt

if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
