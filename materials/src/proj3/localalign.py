from align import align, AffineScoring
from Bio.SubsMat.MatrixInfo import blosum62
from Bio.SeqIO import parse
import sys

def main(queryFile, dbFile, pctId, gapOpen, gapExtend):
    query = parse(queryFile, 'fasta').next()
    db = parse(dbFile, 'fasta')
    score = AffineScoring(blosum62, gapOpen, gapExtend)
    fn = lambda target: align(query, target, score, 'local')
    
    for target in db:
        aln = fn(target)
        if aln.propId * 100.0 > pctId:
            print '>' + target.id
            print '\n'
            print aln
            print '\n\n'

if __name__ == '__main__':
    queryFile = sys.argv[1]
    dbFile = sys.argv[2]
    pctId = int(sys.argv[3])
    gapOpen = int(sys.argv[4])
    gapExtend = int(sys.argv[5])

    main(queryFile, dbFile, pctId, gapOpen, gapExtend)
    
    
