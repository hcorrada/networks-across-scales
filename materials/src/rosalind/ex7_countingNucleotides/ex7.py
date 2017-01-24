import sys

filename = sys.argv[1]
with open(filename, 'r') as f:
    s = f.read().strip()
    tab = dict(A=0, C=0, G=0, T=0)
    for c in s:
        tab[c] = tab[c] + 1
    print " ".join(str(tab[c]) for c in ["A", "C", "G", "T"])

