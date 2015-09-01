# this module provides access to command line arguments
import sys

# get input filename from first command line argument (after script name)
# e.g. $ python pythagoras.py input.txt
filename = sys.argv[1]

# the 'with' construct is useful to deal
# with closing resources etc. (e.g., close the file after opening it)
with open(filename, 'r') as f:
    # read all contents from file (in this case there is only one line)
    # and then 'strip' all leading and trailing whitespace (you should
    # get used to doing this for every line of input you read)
    str = f.read().strip()

    # 'split' the string 'str' into words (any whitespace)
    # then for each word, convert to 'int'
    # assign to variables 'a' and 'b'
    a,b = [int(x) for x in str.split()]

    # compute hypothenuse and print to stdout
    print a**2 + b**2
