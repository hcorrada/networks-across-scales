import sys

def isCG(x):
	return x in ['C', 'G']

def isValid(x):
	return x in ['A', 'C', 'T', 'G']

def calcGC(file):
	totalCount = 0.
	cgCount = 0.
	seqId = None

	for line in file:
		line = line.strip()
		if len(line) == 0:
			continue

		if line[0] == '>':
			if seqId is not None:
				gcContent = 100. * cgCount / totalCount
				print "%s\t%.2f" % (seqId, gcContent)

			seqId = line[1:].split()[0]
			totalCount = 0.
			cgCount = 0.
			continue

		totalCount += sum([isValid(letter) for letter in line])
		cgCount += sum([isCG(letter) for letter in line])

	if totalCount == 0:
		sys.stderr.write('No valid characters in sequence\n')

	elif seqId is not None:
		gcContent = 100. * cgCount / totalCount
		sys.stdout.write("%s\t%.2f\n" % (seqId, gcContent))
	return 0

def main(filename):
	with open(filename, 'r') as file:
		calcGC(file)
	return 0

if __name__ == '__main__':
	if len(sys.argv) < 2:
		print 'not enough arguments'
		sys.exit(1)
	sys.exit(main(sys.argv[1]))

