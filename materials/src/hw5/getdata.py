from Bio import SeqIO, Entrez
Entrez.email = "hcorrada@gmail.com"
id1 = "FN434445.1"
id2 = "FN434454.1"
id3 = "EF566040.1"

handle = Entrez.efetch(db="nucleotide", id=[id1,id2,id3], rettype="gb", retmode="text")

recs = [rec for rec in SeqIO.parse(handle, "genbank")]
SeqIO.write(recs,'data.fa','fasta')
