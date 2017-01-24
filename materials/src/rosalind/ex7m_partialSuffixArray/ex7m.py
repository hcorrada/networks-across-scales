import sys

# construct suffix array from given text
def make_suffix_array(text):
    # setup list of indices
    sa = [i for i in xrange(len(text))]

    # sort the indices by corresponding suffixes
    sa = sorted(sa, key=lambda i: text[i:])
    return sa

def readdat(filename):
    with open(filename, 'r') as f:
        text = f.readline().strip()
        k = int(f.readline().strip())
        return text, k 

def main(filename):
    text, k = readdat(filename)
    sa = make_suffix_array(text)
    sa = zip([i for i in xrange(len(sa))], sa)
    sa = filter(lambda x: x[1] % k == 0, sa)
    for item in sa:
        print "%d,%d" % item
    
if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
    
