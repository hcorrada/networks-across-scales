# Approximate Matching with BWT

Here we implement an approximate matching strategy using
(a) BWT, we use a Burrows-Wheeler transform to perform exact matching, and
(b) use seed-and-check strategy for approximate matching

Usage:

```python
import approximate_matcher

am = ApproximateMatcher(target)
indices = am.get_matches(pattern, d)
```

File `rosalind_approximate_match.py` uses this interface to solve the
Rosalind problem on approximate matching

## Design

`ApproximateMatcher`: computes bwt and suffix array (to make this most space efficient
we would use a partial suffix array and checkpoints), but here we are not doing that. Instead,
you will store a full suffix array, however use precomputed auxiliary data structures to make
BWT matching more efficient.

`get_matches`: argument `d` specifies minimum number of mismatches to allow in approximate
matches. The strategy here will be to use the seed-and-check strategy (as seen in midterm 2).
If `d>0`, then query `pattern` is split into non-overlapping seed kmers as required to ensure that
no approximate matches with d or fewer mismatches are missed when seed kmers are matched exactly.
The bwt is used to find candidate approximate matches by matching each seed kmer exactly. Finally,
candidate approximate matches are checked to ensure that they are indeed within edit distance `d`.

## Things you need to implement

### Seed and check

There is a skeleton of the seed and check strategy in class `SeedChecker` in file `approximate_matcher/seed_and_check/__init__.py`.
You need to implement methods `_make_seeds` and `_within_edit_distance`.

### BWT

There is a skeleton class `BWT` in file `approximate_matcher/bwt/__init__.py`. You need to implement the following
methods:

`_construct(text)`: given a target string construct the Burrows-Wheeler transform and suffix array of string.

`BWT._get_matching_rows(self, pattern)`: returns top and bottom pointers
for rows of the sorted rotations table that contain exact matches to query
string `pattern`. You should use functions `BWT.first_occurence` and
`BWT.count` to implement this. See below. You need to implement this method completely.

`BWT.get_matches(self, pattern)`: returns the indices for positions in target
string that match query string exactly. You need to finish the implementation
of this method.

You also need to implement two functions that pre-process a BWT to create
auxiliary data structures to make BWT exact matching more efficient. These
are in file `approximate_matcher/bwt/preprocess_bwt.py`. Here you need to
complete the implementation of two functions:

`_get_first_occurence_fn(bwt)`: returns a function `f` such that `f(symbol)` is
the index of the first location of `symbol` in first column of rotation
matrix for BWT of target string.

`get_count_fn(bwt)`: returns a function `f` such that `f(symbol, position)`
gives the number of occurrences of `symbol` up to `position` of the BWT of
target string.
