from proj2Lib.ExactMatcher import ExactMatcherNode as ExactMatcherNode
from proj2Lib.ExactMatcher import ExactMatcher as ExactMatcher

# match in O(|T| x |P|) time
class NaiveMatcher(ExactMatcher):
	def preprocessPattern(self, pattern):
		# start at the root node
		currentNode = self.root

		# if we mismatch at the root, advance one position on the target
		currentNode.targetShift = 1

		# add a failure link to the root itself
		currentNode.setFailureLink(self.root)

		for c in pattern:
			# make a new node
			newNode = ExactMatcherNode()

			# add a failure link to the root
			newNode.setFailureLink(self.root)

			# if we mismatch at the new node, we matched the current node
			# so shift at the target one more position to the left
			newNode.targetShift = currentNode.targetShift - 1

			# add the new node as a child indexed by character
			currentNode.setChild(c, newNode)

			# make the new node the current node
			currentNode = newNode

	

