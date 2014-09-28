import sys

filename = sys.argv[1]
with open(filename, 'r') as f:
    s = f.read().strip()
    a,b = map(int, s.split())
    acc = 0
    x = a
    while x <= b:
        if x % 2 == 1:
            acc = acc + x 
        x = x + 1
    print acc
