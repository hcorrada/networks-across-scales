import sys


# make a comparison function based on suffixes of text
def make_cmp_fn(text):
    # comparison function, takes two suffix indices
    def mycmp(x,y):
        # how many positions to check (length of text - largest index)
        n = len(text) - max(x,y)

        # start at the first position of each suffix
        d = cmp(text[x], text[y])
        i = 1

        # keep going while suffixes are equal
        while d == 0 and i < n:
            # compare characters at i'th position
            d = cmp(text[x+i], text[y+i])
            i += 1
        # return the last comparison
        return d
    # return the comparison function
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
    
