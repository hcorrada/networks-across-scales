import sys

def isCG(x):
	return x in ['C', 'G']

def isValid(x):
	return x in ['A', 'C', 'G', 'T']

def getIndex(cgCount, windowSize, nBins):
	index = int(1. * nBins * cgCount / windowSize)
	return index - 1 if index == nBins else index

def outputRes(seqId, bins):
	sys.stdout.write("%s\t%s\n" % (seqId, ','.join(["%d" % x for x in bins])))

def calcGCWindow(file, windowSize=100, nBins=10):
	bins = [0. for i in xrange(nBins)]
	totalCount = 0
	cgCount = 0
	seqId = None
	windowSeq = ['' for i in xrange(windowSize)]

	for line in file:
		line = line.strip()
		if len(line) == 0:
			continue

		if line[0] == '>':
			if seqId is not None:
				outputRes(seqId, bins)

			seqId = line[1:].split()[0]
			bins = [0 for i in xrange(nBins)]
			totalCount = 0
			cgCount = 0
			windowSeq = ['' for i in xrange(windowSize)]
			continue

		if (totalCount  < windowSize):
			need = min(len(line), windowSize - totalCount)
			for letter in line[:need]:
				if not isValid(letter):
					continue

				windowSeq[totalCount] = letter
				totalCount += 1
				cgCount += 1 if isCG(letter) else 0

			if totalCount == windowSize:
				index = getIndex(cgCount, windowSize, nBins)
				bins[index] += 1
				line = line[need:]

		if (totalCount < windowSize):
			continue
		
		nextIndex = 0
		for letter in line:
			otherLetter = windowSeq[nextIndex]
			if not isValid(otherLetter):
				continue

			if isCG(letter) != isCG(otherLetter):
				cgCount += 1 if not isCG(otherLetter) else -1

			index = getIndex(cgCount, windowSize, nBins)
			bins[index] += 1
			windowSeq[nextIndex] = letter
			nextIndex = (nextIndex + 1) % windowSize


	if totalCount < windowSize:
		sys.stderr.write('Not enough valid characters in sequence\n')

	elif seqId is not None:
		outputRes(seqId, bins)
	return 0

def main(filename, windowSize=100, nBins=10):
	with open(filename, 'r') as file:
		calcGCWindow(file, windowSize, nBins)
	return 0

if __name__ == '__main__':
	if len(sys.argv) < 3:
		print 'not enough arguments'
		sys.exit(1)
	sys.exit(main(sys.argv[1],int(sys.argv[2]),int(sys.argv[3])))

