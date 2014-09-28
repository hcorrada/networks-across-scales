import sys

filename = sys.argv[1]
with open(filename, 'r') as f:
    s = f.readline().strip()
    nums = f.readline()

    a,b,c,d = [int(x) for x in nums.split()]
    print s[a:(b+1)] + " " + s[c:(d+1)]

    
