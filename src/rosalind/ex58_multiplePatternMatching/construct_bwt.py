def get_bwt(text):
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

