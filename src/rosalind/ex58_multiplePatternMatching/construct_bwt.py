from bwt_match import preprocess_bwt

def _get_bwt(text):
        # rotate text so that last character of rotation is on the appropriate position
        # i.e. the ith position has the last character of the rotation that
        # starts at the ith position in text (this is not the same order that
        # the rotations are listed on pg. 326 of Compeau & Pevzner
        rotated_text = text[-1] + text[:-1]

        # make tuples of character in bwt and the rotation it corresponds to
        indexed_chars = zip(list(rotated_text), [i for i in xrange(len(text))])
    
        # sort tuples by corresponding rotations
        # the key function return the rotated text starting at appropriate index
        sorted_chars = sorted(indexed_chars, key=lambda x: text[x[1]:] + text[:x[1]])

        # extract the bwt characters after sorting rotations
        bwt_chars, indices = zip(*sorted_chars)

        # turn into a string and return
        return ''.join(list(bwt_chars))


class BWT:
    def __init__(self, text):
        self._bwt = _get_bwt(text)
        self._first_occurence, self._count = preprocess_bwt(self._bwt)

    def move_pointer(self, symbol, index):
        return self._first_occurence[symbol] + self._count(symbol, index)
        
    def move_back(self, index):
        symbol = self._bwt[index]
        return self.move_pointer(symbol, index+1) - 1

    # counts the number of times pattern occurs in text using
    # first and last columns of M(text) matrix reconstructed
    # using preprocess_bwt
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
    
                top = self.move_pointer(symbol, top)
                bottom = self.move_pointer(symbol, bottom + 1) - 1
            else:
                return (top,bottom)
        return (-1,-1)


