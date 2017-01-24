import sys
from Bio import ExPASy
from Bio import SwissProt
import re

def get_processes(id):
    pat = re.compile(r'^P:(?P<process>.*)')
    handle = ExPASy.get_sprot_raw(id)
    record = SwissProt.read(handle)
    entries = filter(lambda x: x[0] == "GO" and pat.match(x[2]) is not None, record.cross_references)
    return [pat.search(entry[2]).group('process') for entry in entries]

filename = sys.argv[1]
with open(filename, 'r') as f:
    id = f.read().strip()
    print "\n".join(get_processes(id))
