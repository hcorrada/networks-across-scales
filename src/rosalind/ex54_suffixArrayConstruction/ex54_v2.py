import sys


def make_cmp_fn(text):
    def mycmp(x,y):
        n = len(text) - max(x,y)
        d = cmp(text[x], text[y])
        i = 1
        
        while d == 0 and i < n:
            d = cmp(text[x+i], text[y+i])
            i += 1
        return d
    return mycmp

# construct suffix array from given text
def make_suffix_array(text):
   
    # setup list of indices
    sa = range(len(text))

    # sort the indices by corresponding suffixes
    sa = sorted(sa, cmp=make_cmp_fn(text))
    return sa

def readdat(filename):
    with open(filename, 'r') as f:
        return f.readline().strip()

def main(filename):
    text = readdat(filename)
    sa = make_suffix_array(text)
    print ", ".join(map(str,sa))
    
if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
    
