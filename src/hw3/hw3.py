from Bio import Entrez
from Bio import SeqIO

Entrez.email = "hcorrada@gmail.com"

handle = Entrez.esearch(db="nuccore",term="257197[BioProject]")
record = Entrez.read(handle)
handle.close()

idlist=record['IdList']

handle = Entrez.efetch(db="nuccore", id=idlist, rettype="gb",retmode="text")
records = SeqIO.parse(handle, "genbank")

with open("ebola.fa", 'w') as f:
    SeqIO.write(records, f, "fasta")
handle.close()
    
rec = SeqIO.read("gbKM233113.txt", format="fasta")
seq = rec.seq

import random

def kdcomposition(text, k, d):
    total_len = 2*k+d
    n = len(text)
    read_pairs = []
    for i in xrange(n-total_len+1):
        kmer1=text[slice(i,i+k)]
        kmer2=text[slice(i+k+d,i+total_len)]
        read_pairs.append(kmer1 + "|" + kmer2)
    random.shuffle(read_pairs)
    return read_pairs

read_pairs = kdcomposition(seq.tostring(), 101, 2000)

# run with pairs

recs=SeqIO.parse("ebola.fa", "fasta")
[(rec.id, rec.seq.tostring() == text) for rec in recs]

read_pairs=kdcomposition(seq.tostring(), 101,20)


read_pairs_20 = kdcomposition(seq.tostring(), 101, 20)
read_pairs_200 = kdcomposition(seq.tostring(), 101, 200)

with open('mystery_20.txt', 'w') as f:
    f.write("\n".join(read_pairs_20))

with open('mystery_200.txt', 'w') as f:
    f.write("\n".join(read_pairs_200))

kmers = kmer_composition(20, seq.tostring())

with open('mystery_single.txt', 'w') as f:
    f.write("\n".join(kmers))
    
