import sys
import re

class gcCounter:
    def __init__(self, s):
        self.theid = s
        self.gc = 0
        self.n = 0
        
    def append(self, s):
        self.n = self.n + len(s)
        self.gc = self.gc + len(re.sub('[AT]', '', s))

    def getId(self): return self.theid
    def getGC(self):
        return float(self.gc) / self.n

    
filename = sys.argv[1]

with open(filename, 'r') as f:
    tab = dict()
    counter = None
    
    while True:
        s = f.readline().strip()
        
        if s is None or len(s) < 1:
            break
        
        if s[0] == ">":
            if counter is not None:
                tab[counter.getId()] = counter.getGC()
            counter = gcCounter(s[1:])
            continue
        counter.append(s)

    if counter is not None:
        tab[counter.getId()] = counter.getGC()
        
    bestid = max(tab, key=lambda x: tab[x])

    print bestid
    print tab[bestid] * 100.
    
