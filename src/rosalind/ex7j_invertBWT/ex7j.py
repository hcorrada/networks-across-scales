import sys
from collections import defaultdict

def tag_occurences(charlist):
    counter = defaultdict(int)
    out = []
    for c in charlist:
        counter[c] += 1
        out.append((c,counter[c]))
    return out

def invert_bwt(bwt):
    bwt = list(bwt)
    last_column = tag_occurences(bwt)    
    first_column = tag_occurences(sorted(bwt))
    
    out = ['$']
    currow = last_column.index(('$',1))
    
    while len(out) < len(bwt):
        entry = first_column[currow]
        out.append(entry[0])
        currow = last_column.index(entry)
    out = ''.join(out)
    return out[1:] + out[0]

def readdat(filename):
    with open(filename, 'r') as f:
        return f.readline().strip()

def main(filename):
    bwt = readdat(filename)
    text = invert_bwt(bwt)
    print text

if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
