import sys

filename = sys.argv[1]
with open(filename, 'r') as f:
    s = f.read().strip()
    out = s.replace('T', 'U')
    print out
