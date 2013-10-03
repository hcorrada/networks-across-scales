class ExactMatcherNode(object):
	def __init__(self, root=False):
		self.children = dict() # children nodes labeled by character
		self.depth = 0 # node depth in graph

		self.failureLink = None # node to go to on mismatch

		self.targetShift = 0 # how to update target pointer after mismatch

		self.numQueries = 0 # number of matching queries against node

	def printNode(self, full=True):
		childrenChars = ','.join(self.children.keys())
		if not full:
			return '<ExactMatcherNode: %s, depth %d>' % (childrenChars, self.depth)

		pad = ''.join(['|' for i in xrange(self.depth)])
		failLink = self.failureLink.printNode(full=False)
		return '%s<ExactMatcherNode:  %s, failLink: %s, targetShift %d, depth %d>\n' % (pad, childrenChars, failLink, self.targetShift, self.depth)

	def __repr__(self):
		return self.printNode(full=True)

	def setFailureLink(self, node):
		self.failureLink = node

	def setChild(self, c, node):
		if self.isMatch(c):
			return False

		node.depth = self.depth + 1
		self.children[c] = node
		return True

	def getChildLabels(self):
		return self.children.keys()

	def getChild(self, index):
		if self.isLeaf():
			return None
		c = self.getChildLabels()[0]
		return self.children[c]

	def isLeaf(self): return len(self.children) == 0
	def isRoot(self): return self.depth == 0

	def isMatch(self, c):
		self.numQueries += 1
		return c in self.children

	def getTransition(self, c):
		if not self.isMatch(c):
			# mismatch or leaf, treated equally: follow failure link and shift target accordingly
			return (self.failureLink, self.targetShift)

		# match, follow transition and shift target by one
		return (self.children[c], 1)

		
class ExactMatcher(object):
	def __init__(self, pattern, reportFormat = 'Found match at position {0}'):
		assert len(pattern) > 0
		self.root = ExactMatcherNode(root=True)
		self.patternLength = len(pattern)
		self.preprocessPattern(pattern)
		self.reportFormat = reportFormat

	def reportMatch(self, i):
		print self.reportFormat.format(i - self.getPatternLength())

	def getPatternLength(self): return self.patternLength

	def matchTarget(self, target):
		assert len(target) >= self.getPatternLength()

		# start at the root
		currentNode = self.root

		# start in the first position of the target
		i = 0  
		while i < len(target):
			c = target[i]
#			print 'i:', i, c

#			print 'curNode:', currentNode

			# get the current target character

			# check if this is a leaf (pattern is exhausted)
			if currentNode.isLeaf():
				self.reportMatch(i)

			# try to matching it to an exiting edge on the current node
			# if current node is a leaf, this function updates correctly
			(currentNode, targetShift) = currentNode.getTransition(c)

#			print 'shift', targetShift, '\n'
			i += targetShift

	# stat reporting
	def walker(self, fun):
		def walkerHelper(node):
			out = fun(node)

			if node.isLeaf():
				return out

			for childNode in node.children.itervalues():
				out += walkerHelper(childNode)
			return out

		return walkerHelper(self.root)

	def getNumQueries(self):
		return self.walker(lambda node: node.numQueries)

	def getTotalTime(self):	return self.getNumQueries()

	def __repr__(self):
		return self.walker(lambda node: repr(node))



