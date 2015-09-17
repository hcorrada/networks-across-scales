import sys

from Bio import Entrez
from Bio import SeqIO

# download fasta records for given ids
# input:
#   ids: a list of ids <strings>
# output:
#   a list of records as SeqIO objects
def download_fasta_records(ids):
    Entrez.email = "hcorrada@gmail.com"
    id = ','.join(ids)
    handle = Entrez.efetch(db="nucleotide", id=ids, rettype="fasta")
    records = list(SeqIO.parse(handle, "fasta"))
    return records

def find_shortest_sequence(records):
    return sorted(records, key=lambda x: len(x.seq))[0]

def print_fasta_record(rec):
    SeqIO.write(rec, sys.stdout, "fasta")

filename = sys.argv[1]
with open(filename, 'r') as f:
    ids = f.readline().strip().split()
    records = download_fasta_records(ids)
    rec = find_shortest_sequence(records)
    print_fasta_record(rec)
