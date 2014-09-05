import sys

filename = sys.argv[1]
with open(filename, 'r') as f:
    s = f.read().strip().split()
    tab = dict()
    for word in s:
        if word not in tab:
            tab[word] = 1
        else:
            tab[word] = tab[word] + 1
    for word, cnt in tab.iteritems():
        print word, cnt

    
