from Bio import Entrez, SeqIO

Entrez.email = "hcorrada@gmail.com"
handle = Entrez.efetch(db="nucleotide", id="AL111168.1", rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()

SeqIO.write(record, "campy.fa", "fasta")


record = SeqIO.read("campy.fa", "fasta")

def skew_vec(genome):
    res = [0]
    skew = 0
    for i in xrange(len(genome)):
        c = genome[i]

        if c == 'C':
            skew = skew - 1
        elif c == 'G':
            skew = skew + 1
        res.append(skew)
    return res

vec = skew_vec(record)

