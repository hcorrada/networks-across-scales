# A node in an exact matching automata
class ExactMatcherNode(object):
	def __init__(self, root=False):
		self.children = dict() # children nodes labeled by character
		self.depth = 0 # node depth in graph

		self.failureLink = None # node to go to on mismatch

		self.targetShift = 0 # how to update target pointer after mismatch

		self.numQueries = 0 # number of matching queries against node

		self.outputLink = None

	# pretty print node information
	def printNode(self, full=True):
		childrenChars = ','.join(self.children.keys())
		if not full:
			return '<ExactMatcherNode: %s, depth %d>' % (childrenChars, self.depth)

		pad = ''.join(['|' for i in xrange(self.depth)])
		failLink = self.failureLink.printNode(full=False)
		return '%s<ExactMatcherNode:  %s, failLink: %s, targetShift %d, depth %d>\n' % (pad, childrenChars, failLink, self.targetShift, self.depth)

	def __repr__(self):
		return self.printNode(full=True)

	# set the failure link
	def setFailureLink(self, node):
		self.failureLink = node

	# set the output link
	def setOutputLink(self, node):
		self.outputLink = node
		
	# set child on outgoing edge labeled with character c
	def setChild(self, c, node):
		if self.isMatch(c):
			return False

		node.depth = self.depth + 1
		self.children[c] = node
		return True

	# get labels of exiting edges
	def getChildLabels(self):
		return self.children.keys()

	# get child at index (integer) of outgoing edge set
	def getChild(self, index):
		if self.isLeaf():
			return None
		c = self.getChildLabels()[0]
		return self.children[c]

	# is this a leaf?
	def isLeaf(self): return len(self.children) == 0

	# is this the root?
	def isRoot(self): return self.depth == 0

	# match outgoing edge to character c
	# THIS IS THE FUNCTION THAT KEEP TRACKS OF CHARACTER COMPARISONS
	# MAKE SURE THIS IS USED WHENEVER A COMPARISON IS MADE
	def isMatch(self, c):
		self.numQueries += 1
		return c in self.children

	# Find out how to transition out of this node when matching character c
	# returns the next node to visit and how to shift a pointer into the target string
	def getTransition(self, c):
		if not self.isMatch(c):
			# mismatch or leaf, treated equally: follow failure link and shift target accordingly
			return (self.failureLink, self.targetShift)

		# match, follow transition and shift target by one
		return (self.children[c], 1)

		
# abstract class for exact matcher
class ExactMatcher(object):
	def __init__(self, pattern, reportFormat = 'Found match at position {0}'):
		assert len(pattern) > 0
		self.root = ExactMatcherNode(root=True) # the root node
		self.patternLength = len(pattern) # the length of the pattern
		self.preprocessPattern(pattern) # preprocess the pattern, subclasses need to implement this
		self.reportFormat = reportFormat # format string to use when reporting a match


	# use the given format string to report a match
	# i is the position on the target where matching pattern ends
	def reportMatch(self, i):
		print self.reportFormat.format(i - self.getPatternLength())

	# get pattern length
	def getPatternLength(self): return self.patternLength

	# generic implementation of matching through automata
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

	# helper function that walks the tree to obtains information about automata
	def walker(self, fun):
		def walkerHelper(node):
			out = fun(node)

			if node.isLeaf():
				return out

			for childNode in node.children.itervalues():
				out += walkerHelper(childNode)
			return out

		return walkerHelper(self.root)

	# get the total number of character comparisons
	def getNumQueries(self):
		return self.walker(lambda node: node.numQueries)

	# same as above
	def getTotalTime(self):	return self.getNumQueries()

	# pretty print the automata
	def __repr__(self):
		return self.walker(lambda node: repr(node))



