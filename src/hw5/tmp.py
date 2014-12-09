from Bio import SeqIO
from kmer_composition import kmer_composition
from bwt import BWT
from suffix_array import SuffixArray
from approximate_match import find_approximate_matches
from pileup import PileUp
from simreads import simulate_reads

with open('data.fa', 'r') as f:
    recs = [str(rec.seq) for rec in SeqIO.parse(f,'fasta')]

reference = recs[0]
reads = simulate_reads(recs[1], 50, 10, .05)

def make_pileup(kmers, reference):
    p = PileUp(reference)
    bwt = BWT(reference + '$')
    sa = SuffixArray(reference + '$')

    for kmer in kmers:
        indices = find_approximate_matches(kmer, bwt, sa, 3, reference)
        if len(indices) == 0:
            continue
        
        indices = list(zip(*indices)[0])
        p.insert(indices, kmer)
    return p

p = make_pileup(reads, reference)
p.print_mismatches(.5)





    
