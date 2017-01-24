from preprocess_bwt import get_first_occurence, count_occurences, count_occurences_checkpoints

# comparison function for construction
# of burrows-wheeler transform
def make_cmp(text):
    n = len(text)
    
    def _cmp(x,y):
        d = cmp(text[x], text[y])
        i = x
        j = y
        k = 0
        while d == 0 and k < n:
            i = (i + 1) % n
            j = (j + 1) % n

            d = cmp(text[i], text[j])
        return d
    return _cmp

# compute the bwt of given text
def _get_bwt(text):
        # setup rotation index array
    indices = range(len(text))

    # sort index array by string rotations
    # uses comparison function defined above
    indices = sorted(indices, cmp=make_cmp(text))

    # get the last character of each rotation
    bwt = [text[i-1] for i in indices]

    # turn into a string and return
    return ''.join(bwt)

# class encapsulating bwt operations
class BWT:
    # initialize object use a checkpoint table if checkpoints > 0
    def __init__(self, text, checkpoints=0):
        # compute bwt of given text
        self._bwt = list(_get_bwt(text))

        # find first occurences in first column of rotation matrix
        self._first_occurence = get_first_occurence(self._bwt)

        # get the count function
        if checkpoints > 0:
            # using checkpoints
            self._count = count_occurences_checkpoints(self._bwt, checkpoints)
        else:
            # without checkpoints
            self._count = count_occurences(self._bwt)

    # find pointer in first column of last occurence of symbol
    # before 'index' in last column
    def move_pointer(self, symbol, index):
        return self._first_occurence[symbol] + self._count(symbol, index)

    # find pointer in first column of suffix starting with
    # symbol at index in last column
    def move_back(self, index):
        symbol = self._bwt[index]
        return self.move_pointer(symbol, index+1) - 1

    # find range of matching rows in rotation matrix
    # containing pattern
    # returns (-1,-1) if no matches
    def find_matches(self, pattern):
        # have top and bottom pointers at first and last
        # row of M matrix
        top = 0
        bottom = len(self._bwt) - 1

        # while there are rows to be checked
        while top <= bottom:
            # and patterns symbols to be checked
            if len(pattern) > 0: 
                symbol = pattern[-1] # check the last symbol in pattern
                pattern = pattern[:-1] # remove last symbol in pattern

                # update top and bottom pointers
                top = self.move_pointer(symbol, top)
                bottom = self.move_pointer(symbol, bottom + 1) - 1
            else:
                return (top,bottom)
        return (-1,-1)


