import sys

def getGC(s):
    cnt = 0
    for c in s:
        if c == "G" or c == "C":
            cnt = cnt + 1
    return cnt, len(s)

filename = sys.argv[1]

bestId = ''
maxGC = 0.

with open(filename, 'r') as f:
    gc = 0
    n = 0
    strid = ''
    
    while True:
        s = f.readline().strip()
        
        if s is None or len(s) < 1:
            break
        
        if s[0] == ">":
            if n > 0:
                gcContent = float(gc) / n
                print strid, gcContent
                if gcContent > maxGC:
                    bestId = strid
                    maxGC = gcContent
            strid = s[1:]
            gc = 0
            n = 0
        else:
            curGC, curN = getGC(s)
            gc = gc + curGC
            n = n + curN
    if n > 0:
        gcContent = float(gc) / n
        print strid, gcContent
        if gcContent > maxGC:
            bestId = strid
            maxGC = gcContent
        
    print bestId
    print maxGC * 100.
