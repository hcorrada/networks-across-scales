import sys

class Path:
    def __init__(self):
        self.__slot__ = None

    def get_string(self):
        return None

def build_pairedend_path(paired_kmers, k, d):
    return Path()

def readdat(filename):
    return None, None, None

def main(filename):
    paired_kmers, k, d = readdat(filename)
    path = build_pairedend_path(paired_kmers, k, d)
    string = path.get_string()
    print string

# this is here so this plays nicely with ipython %loadpy magic
if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
