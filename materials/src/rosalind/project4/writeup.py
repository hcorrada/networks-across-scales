from Bio import SeqIO, Seq
from approximate_matcher import ApproximateMatcher
from pileup import PileUp, PileUpStrategy
import itertools

d = 3
with open('data/reference.fa', 'r') as reference_file, open('data/reads.fa', 'r') as reads_file:
    reference = SeqIO.read(reference_file, 'fasta')
    print reference

    amatcher = ApproximateMatcher(str(reference.seq))
    pileup = PileUp(reference)

    reads = list(SeqIO.parse(reads_file, 'fasta'))

    for read in reads:
        pattern = read.seq
        indices = amatcher.get_matches(pattern, d)

        if len(indices) > 0:
            pileup.insert(indices, pattern)

    pileup.print_mismatches(0.5)

    # position 822 C->T
    print reference.seq.translate()[270:280]
    consensus = str(reference.seq)
    consensus = consensus[:822] + 'T' + consensus[823:]
    consensus = Seq.Seq(consensus)

    print consensus.translate()[270:280]
