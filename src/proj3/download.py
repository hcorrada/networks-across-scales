from Bio import Entrez
Entrez.email = "hcorrada@umiacs.umd.edu"

ids=['ADJ40477.1','AGM53846.1','AAA91327.1','ABQ10080.1','ABQ09751.1','ADI99553.1','ADI99633.1']

handle = Entrez.efetch(db="protein", id=ids,rettype="gb",retmode="text")

from Bio import SeqIO
records = SeqIO.parse(handle,"gb")

fp = open('flustrains.fa', 'w')
SeqIO.write(records,fp,'fasta')
fp.close()

records = SeqIO.parse('flustrains.fa', 'fasta')
query = records.next()

SeqIO.write(query, 'flustrains-query.fa', 'fasta')
SeqIO.write(records, 'flustrains-db.fa','fasta')

from align import align, AffineScoring
x=SeqIO.parse('flustrains-query.fa','fasta').next()
y=SeqIO.parse('flustrains-db.fa','fasta').next()

from SeqIO.SubsMat.MatrixInfo import blosum62
ascore=AffineScoring(blosum62,-11,-1)

res=align(x,y,ascore,'local')

query=SeqIO.parse('flustrains-query.fa','fasta').next()
db=SeqIO.parse('flustrains-db.fa','fasta')
print "Query {0.description}".format(query)
for target in db:
    aln=align(query,target,ascore,'local')
    print "Target {0.description}: pctId {1.propId:.1%}, Score {1.score}".format(target,aln)


nejmids = ['FN434454','FN434445']
nejmids2=['CBG91923.1','CBG91982.1']
handle = Entrez.efetch(db="protein",id=nejmids2,rettype="gb",retmode="text")

records=SeqIO.parse(handle,"gb")
SeqIO.write(records,'flustrains-nejm.fa','fasta')

records=SeqIO.parse('flustrains-nejm.fa','fasta')
for target in records:
    print align(query,target,ascore,'local')

import re
from Bio import SeqIO

pat = re.compile('.*\[Influenza A virus \((.*)\)\]')
with open('flustrains-db.fa','r') as handle:
    for rec in SeqIO.parse(handle,'fasta'):
        match = pat.match(rec.description)
        print match.group(1)
        
