import sys
from Bio import Entrez

def count_genbank_entries(genus, start_date, end_date):
    Entrez.email = "hcorrada@gmail.com"
    handle = Entrez.esearch(db="nucleotide", term="%s[Organism" % genus, datetype="pdat", mindate=start_date, maxdate=end_date)
    record = Entrez.read(handle)
    return record['Count']

filename = sys.argv[1]
with open(filename, 'r') as f:
    genus = f.readline().strip()
    start_date = f.readline().strip()
    end_date = f.readline().strip()
    print count_genbank_entries(genus, start_date, end_date)
