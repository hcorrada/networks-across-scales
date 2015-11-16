import sys

def construct_suffix_array(text):
    return range(len(text))

def find_matches(pattern, sa):
    return range(3)

def readdat(filename):
    return "abcd", ["ab", "cd"]

def main(filename):
    text, patterns = readdat(filename)
    sa = construct_suffix_array(text)

    res = []
    for pattern in patterns:
        res += find_matches(pattern, sa)
    print " ".join(map(str, res))

if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
