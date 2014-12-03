import sys

def get_bwt(text):
    # setup rotation index array
    indices = range(len(text))

    # sort index array by string rotations
    indices = sorted(indices, key=lambda i: text[i:] + text[:i])

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
