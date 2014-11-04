import sys
import numpy as np

def dpchange(money, coins):
    minnumcoins = np.zeros(money+1)
    minnumcoins[0] = 0
    for m in xrange(1,money+1):
        minnumcoins[m] = float('inf')
        for coin in coins:
            if m >= coin:
                tmp = minnumcoins[m-coin] + 1
                if tmp < minnumcoins[m]:
                    minnumcoins[m] = tmp
    return int(minnumcoins[money])

def readdat(filename):
    with open(filename, 'r') as f:
        money = int(f.readline().strip())
        coins = map(int, f.readline().strip().split(','))
    return money, coins

def main(filename):
    money, coins = readdat(filename)
    out = dpchange(money, coins)
    print out
    
if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
