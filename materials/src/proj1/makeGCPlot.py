from pylab import *

# read files 
files = ['bacteria.out', 'Celegans.out', 'Athaliana.out', 'Hsapiens.out']
ids = []
gc = []

for file in files:
	with open(file, 'r') as f:
		for line in f:
			tmp = line.split('\t')
			ids.append(tmp[0])
			gc.append(float(tmp[1]))

barh(range(len(ids)), gc)
savefig("fig1.pdf", format="pdf")



files = ['bacteria_win.out', 'Celegans_win.out', 'Athaliana_win.out', 'Hsapiens_win.out']
ids = []
gc = []

for file in files:
	with open(file, 'r') as f:
		for line in f:
			tmp = line.split('\t')
			ids.append(tmp[0])
			gc.append([int(x) for x in tmp[1].split(',')])

for i in range(len(ids)):
	h = gc[i]
	figure(i)
	bar(range(len(h)), h)
	xlim(0,len(h))
	title(ids[i])
	savefig("fig2_%s.pdf" % ids[i].strip(), format="pdf")

