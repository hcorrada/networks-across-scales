from Bio import SeqIO

db = SeqIO.parse("bigtest-db.fa","fasta")
out = []

i = 0
for seq in db:
    if i == 10:
        break

    out.append(seq)
    i += 1

outf = open("bigtest2-db.fa",'w')
SeqIO.write(out,outf,'fasta')
outf.close()
