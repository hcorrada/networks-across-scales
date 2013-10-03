# match in O(|T| + |P|) 
from proj2Lib.ExactMatcher import ExactMatcher as ExactMatcher
from proj2Lib.ExactMatcher import ExactMatcherNode as ExactMatcherNode

class KMPMatcher(ExactMatcher):
	def findFailNode(self, node, x):
		currentNode = node.failureLink # v
		while (not currentNode.isMatch(x) and not currentNode.isRoot()):
			currentNode = currentNode.failureLink
		if currentNode.isMatch(x):
			return currentNode.getChild(0)
		else:
			return self.root

	def preprocessPattern(self, pattern):
		currentNode = self.root

		currentNode.targetShift = 1
		currentNode.setFailureLink(self.root)

		for i in xrange(len(pattern)):
			c = pattern[i]
			newNode = ExactMatcherNode()

			failNode = self.findFailNode(currentNode, c)			
			newNode.setFailureLink(failNode)
			newNode.targetShift = 0

			currentNode.setChild(c, newNode)
			currentNode = newNode

		# optimize failure links
		currentNode = self.root
		while not currentNode.isLeaf():
			nextNode = currentNode.getChild(0)
			v = currentNode.failureLink
			x = '$' if nextNode.isLeaf() else nextNode.getChild(0)
			if v.isMatch(x):
				nextNode.setFailureLink(v.failureLink)
			currentNode = nextNode

