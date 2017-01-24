import sys

from Bio import Entrez
from Bio import SeqIO

# download fasta records for given ids
# input:
#   ids: a list of ids <strings>
# output:
#   a list of records as Seq objects
def download_fasta_records(ids):
    Entrez.email = "hcorrada@gmail.com"
    id = ','.join(ids)
    handle = Entrez.efetch(db="nucleotide", id=ids, rettype="fasta")
    records = list(SeqIO.parse(handle, "fasta"))
    return records

# find the shortest sequence within a list of Seq records
# input:
#   records: a list of Seq objects
# output:
#   a Seq for the record with the shortest sequence
def find_shortest_sequence(records):
    return sorted(records, key=lambda x: len(x.seq))[0]

# print a Seq record to stdout as fasta
# input:
#   rec: a Seq record
# output:
#   NONE
def print_fasta_record(rec):
    SeqIO.write(rec, sys.stdout, "fasta")

filename = sys.argv[1]
with open(filename, 'r') as f:
    ids = f.readline().strip().split()
    records = download_fasta_records(ids)
    rec = find_shortest_sequence(records)
    print_fasta_record(rec)
