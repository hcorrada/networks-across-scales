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

filename = sys.argv[1]
with open(filename, 'r') as f:
    pattern = f.read().strip()
    print translate(pattern)
