from Bio import SeqIO
from kmer_composition import kmer_composition
from approximate_match import find_approximate_matches

recs = [rec for rec in SeqIO.parse('data.fa','fasta')]
kmers = kmer_composition(50, str(recs[0].seq))

reference = str(recs[1].seq)
bwt = BWT(reference + '$', checkpoints=10)
sa = PartialSuffixArray(reference + '$', 10, bwt)

matches = {}
for kmer in kmers:
    top, bottom = bwt.find_matches(kmer)
    if top == -1:
        continue
    matches[kmer] = sa.get_match_indices(top, bottom)

import operator
indices = sorted(map(operator.itemgetter(0), matches.values()))

import numpy as np
indices = np.array(indices)
d=indices[1:]-indices[:-1]

