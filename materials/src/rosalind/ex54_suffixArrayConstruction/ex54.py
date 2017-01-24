import sys

# construct suffix array from given text
def make_suffix_array(text):
    # setup list of indices
    sa = range(len(text))

    # sort the indices by corresponding suffixes
    sa = sorted(sa, key=lambda i: text[i:])
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
    
