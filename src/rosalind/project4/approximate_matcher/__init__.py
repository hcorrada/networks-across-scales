from bwt import BWT
from seed_and_check import SeedChecker

# class encapsulating approximate matching with
# seed and check strategy, using BWT for exact matching
# of seeds
#
# Fields:
#   _text: target string, with '$' character appended
#   _bwt: object of class BWT, used for exact matching of seeds
class ApproximateMatcher:
    def __init__(self, target):
        self._text = target + '$'
        self._bwt = BWT(self._text)

    # return indices in target that contain
    # matches of string pattern with up to d
    # mismatches
    def get_matches(self, pattern, d):
        # initialze seed and check object
        seed_checker = SeedChecker(pattern, d)

        # for each seed k-mer in pattern
        for seed, seed_index in seed_checker.enumerate():
            # find exact matches of seed using BWT
            indices = self._bwt.get_matches(seed)
            # add candidate approximate matches based on
            # seed exact matches
            seed_checker.add_candidates(indices, seed_index)
        # verify that candidate approximate matches are within
        # minimum edit distance, and return final matches
        matches = seed_checker.filter_candidates(self._text)
        return matches
