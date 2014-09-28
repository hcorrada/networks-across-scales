import sys

filename = sys.argv[1]
with open(filename, 'r') as f:
    i = 1
    for line in f:
        if i % 2 == 0:
            print line.strip()
        i = i +1

