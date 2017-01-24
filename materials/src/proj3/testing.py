alphabet = ['A','C','G','T']
score = dict([((x,y), 2 if x==y else -1) for x in alphabet for y in alphabet])


from Bio import SeqIO

x = SeqIO.parse("smalltest-query.fa","fasta").next()
y = SeqIO.to_dict(SeqIO.parse("smalltest-db.fa","fasta"))['exact']

#mat = buildLinearMatrix(x,y,score,-2)
#getLinearAlignment(mat,x,y)

from align import LinearScoring, align, Aligner
lscore = LinearScoring(score,-5)

def runit0():
    mat = Aligner(x,y,lscore,'local')
    mat.init()
    mat.fill()
    print mat
    print mat.backtrace()

def runit():
    recs = SeqIO.parse("smalltest-db.fa", "fasta")
    for y in recs:
        print y.id
        print x.seq
        print y.seq + '\n'
        print 'linear global'
        print align(x,y,LinearScoring(score,-5),alnKind='global') 

        print 'linear local'
        print align(x,y,LinearScoring(score,-5),alnKind='local')

        print 'affine global'
        print align(x,y,AffineScoring(score,-5, -1),alnKind='global') 

        print 'affine local'
        print align(x,y,AffineScoring(score,-5, -1),alnKind='local')

        print ''

