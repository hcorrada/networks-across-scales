from path import Path
from string_builder import build_string
from graph import Node
from paired_kmer import PairedKmer

k=4
d=2

path = Path()

path.append(Node(PairedKmer("GACC|GCGC"), 0))
path.append(Node(PairedKmer("ACCG|CGCC"), 1))
path.append(Node(PairedKmer("CCGA|GCCG"), 2))
path.append(Node(PairedKmer("CGAG|CCGG"), 3))
path.append(Node(PairedKmer("GAGC|CGGA"), 4))

print path
string = build_string(path, k, d)
print string
