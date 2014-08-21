from Bio import SeqIO
from Bio.SubsMat.MatrixInfo import blosum62

x = SeqIO.parse("bigtest-query.fa","fasta").next()
y = SeqIO.parse("bigtest-db.fa","fasta").next()

#mat = buildLinearMatrix(x,y,score,-2)
#getLinearAlignment(mat,x,y)

from align import AffineScoring, align, Aligner
ascore = AffineScoring(blosum62,-11, -1)

def runit0():
    mat = Aligner(x,y,ascore,'global')
    mat.init()
    mat.fill()
    return mat.backtrace()

def runit():
    x = SeqIO.parse("smalltest-query2.fa", "fasta").next()
    recs = SeqIO.parse("smalltest-db2.fa", "fasta")
    for y in recs:
        print y.id
        print x.seq
        print y.seq + '\n'
        print 'linear global'
        print align(x,y,LinearScoring(blosum62,-5),alnKind='global') 

        print 'linear local'
        print align(x,y,LinearScoring(blosum62,-5),alnKind='local')

        print 'affine global'
        print align(x,y,AffineScoring(blosum62,-5, -1),alnKind='global') 

        print 'affine local'
        print align(x,y,AffineScoring(blosum62,-5, -1),alnKind='local')

        print ''

