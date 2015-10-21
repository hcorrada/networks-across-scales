import sys
import numpy as np

# make change using dynamic programming
# return the minimum number of coins of given
# denominations that add to money
#
# input:
#   money: integer, total amount of money required
#   coins: list of integers, coin denominations
#
def dp_change(money, coins):
    # initialize 1-dimensional dynamic programming table
    # to infinity in all entries
    min_num_coins = np.inf * np.ones(money+1)

    # the minimum number of coins needed to make 0 money
    min_num_coins[0] = 0

    # now fill in the rest of the table
    for cur_money in xrange(1, money+1):
        # see what happens if we add one of each coin
        for coin in coins:
            # if coin is larger than the amount of
            # current money we need, we cant use it
            if cur_money >= coin:
                # let's add one of these coins to the minimum
                # number of coins needed to make money - coin
                cur_num_coins = min_num_coins[cur_money-coin] + 1

                # if this is better than what we have found so far
                # store it in the table
                if cur_num_coins < min_num_coins[cur_money]:
                    min_num_coins[cur_money] = cur_num_coins
    # now the answer will be in the last entry
    # of the table, cast to int since np.inf is float
    return int(min_num_coins[money])

def readdat(filename):
    with open(filename, 'r') as f:
        money = int(f.readline().strip())
        coins = map(int, f.readline().strip().split(','))
    return money, coins

def main(filename):
    money, coins = readdat(filename)
    out = dp_change(money, coins)
    print out

# ipython %loadpy trick
if __name__ == '__main__' and 'get_ipython' not in dir():
    filename = sys.argv[1]
    main(filename)
