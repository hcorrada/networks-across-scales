import sys

alphabet = ['A', 'C', 'G', 'U']

translate_tab = {}
translate_tab["A"] = dict.fromkeys(alphabet, {})
translate_tab['A']['A'] = dict.fromkeys(['C','U'], 'N')
translate_tab['A']['A'].update(dict.fromkeys(['A','G'], 'K'))

translate_tab['A']['C'] = dict.fromkeys(alphabet, 'T')

translate_tab['A']['G'] = dict.fromkeys(['U','C'], 'S')
translate_tab['A']['G'].update(dict.fromkeys(['A','G'], 'R'))

translate_tab['A']['U'] = dict.fromkeys(['U','C','A'], 'I')
translate_tab['A']['U']['G'] = 'M'

translate_tab['C'] = dict.fromkeys(alphabet, {})
translate_tab['C']['A'] = dict.fromkeys(['U','C'], 'H')
translate_tab['C']['A'].update(dict.fromkeys(['A','G'], 'Q'))

translate_tab['C']['C'] = dict.fromkeys(alphabet, 'P')
translate_tab['C']['G'] = dict.fromkeys(alphabet, 'R')
translate_tab['C']['U'] = dict.fromkeys(alphabet, 'L')

translate_tab['G'] = dict.fromkeys(alphabet, {})
translate_tab['G']['A'] = dict.fromkeys(['U','C'], 'D')
translate_tab['G']['A'].update(dict.fromkeys(['A','G'], 'E'))

translate_tab['G']['C'] = dict.fromkeys(alphabet, 'A')
translate_tab['G']['G'] = dict.fromkeys(alphabet, 'G')
translate_tab['G']['U'] = dict.fromkeys(alphabet, 'V')

translate_tab["U"] = dict.fromkeys(alphabet, {})
translate_tab['U']['A'] = dict.fromkeys(['C','U'], 'Y')
translate_tab['U']['A'].update(dict.fromkeys(['A','G'], '*'))

translate_tab['U']['C'] = dict.fromkeys(alphabet, 'S')

translate_tab['U']['U'] = dict.fromkeys(['U','C'], 'F')
translate_tab['U']['U'].update(dict.fromkeys(['A','G'], 'L'))

translate_tab['U']['G'] = dict.fromkeys(['U','C'], 'C')
translate_tab['U']['G']['G'] = 'W'
translate_tab['U']['G']['A'] = '*'

def translate_codon(codon):
  x,y,z = list(codon)
  res = translate_tab[x][y][z]
  return res if res != '*' else ''


def translate(pattern):
    res = ''
    for i in xrange(0,len(pattern),3):
        res += translate_codon(pattern[i:i+3])
    return res

def transcribe(pattern):
    out = pattern.replace('T','U')
    return out

def find_matches(text, peptide):
    k = len(peptide) * 3
    n = len(text)
    res = []
    
    for i in xrange(3):
        for j in xrange(i, n - k + 1, 3):
            candidate_peptide = translate(text[j:j+k])
            if (candidate_peptide == peptide):
                res.append(j)
    return res

complement = dict(A="T", C="G", G="C", T="A")

def revcomp(text):
    out = list(text)
    out.reverse()
    for i in xrange(len(out)):
        out[i] = complement[out[i]]
    return ''.join(out)

def peptide_encoding(text, peptide):
    rna = transcribe(text)
    res = find_matches(rna, peptide)
    rna = transcribe(revcomp(text))
    rev_comp_res = find_matches(rna, peptide)
    n = len(text)
    k = len(peptide) * 3
    rev_comp_res = [n - i - k for i in rev_comp_res]
    res += rev_comp_res
    return [text[i:i+k] for i in res]

filename = sys.argv[1]
with open(filename, 'r') as f:
    text = f.readline().strip()
    peptide = f.readline().strip()
    res = peptide_encoding(text, peptide)
    print "\n".join(res)
