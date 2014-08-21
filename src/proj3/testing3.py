from Bio import SeqIO
from Bio.SubsMat.MatrixInfo import blosum62
from align import AffineScoring, align, Aligner

query = SeqIO.parse("bigtest-query.fa","fasta").next()
db = SeqIO.parse("bigtest-db.fa","fasta")

score=AffineScoring(blosum62,-11,-1)

fn = lambda target: align(query, target,score,'local')

i=0
for target in db:
    if i == 10:
        break
    
    print target.id, 
    aln = fn(target)
    print aln.propId
    i += 1

    
    

