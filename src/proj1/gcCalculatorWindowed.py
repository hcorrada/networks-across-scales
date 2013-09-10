import sys
from collections import deque
import matplotlib.pyplot as plt
import numpy as np

def isCG(x):
	return x in ['C', 'G']

def getIndex(cgCount, windowSize, nBins):
	index = int(1. * nBins * cgCount / windowSize)
	return index - 1 if index == nBins else index

def getDensity(bins):
	s = float(sum(bins))
	dens = [x/s for x in bins]
	assert(sum(dens) == 1.)
	return dens

def outputRes(seqId, bins):
	print seqId,"\t", ','.join(["%d" % x for x in bins])

def calcGCWindow(file, windowSize=100, nBins=10):
	bins = [0. for i in xrange(nBins)]
	totalCount = 0
	cgCount = 0
	seqId = None
	windowSeq = deque()

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
			windowSeq = deque()
			continue

		if (len(windowSeq) < windowSize):
			need = windowSize - len(windowSeq)
			for letter in line[:need]:
				windowSeq.append(letter)
				cgCount += 1 if isCG(letter) else 0
			index = getIndex(cgCount, windowSize, nBins)
			bins[index] += 1
			line = line[need:]

		if (len(windowSeq) < windowSize):
			continue

		assert(len(windowSeq) == windowSize)
		for letter in line:
			otherLetter = windowSeq.popleft()

			if isCG(letter) != isCG(otherLetter):
				cgCount += 1 if not isCG(otherLetter) else -1

			index = getIndex(cgCount, windowSize, nBins)
			bins[index] += 1
			windowSeq.append(letter)
			

	if seqId is not None:
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

