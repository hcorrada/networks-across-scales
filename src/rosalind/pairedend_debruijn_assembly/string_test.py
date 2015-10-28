from path import Path
from string_builder import build_string
from graph import Node
from paired_kmer import PairedKmer

k=4
d=2

path = Path()

path.append(Node(PairedKmer("GACC|GCGC")))
path.append(Node(PairedKmer("ACCG|CGCC")))
path.append(Node(PairedKmer("CCGA|GCCG")))
path.append(Node(PairedKmer("CGAG|CCGG")))
path.append(Node(PairedKmer("GAGC|CGGA")))

print path
string = build_string(path, k, d)
print string
