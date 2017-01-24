from proj2Lib.ExactMatcher import ExactMatcher as ExactMatcher
from proj2Lib.ExactMatcherNode import ExactMatcherNode as ExactMatcherNode

# UPDATE THIS FILE TO IMPLEMENT THE AHO-CORASICK ALGORITHM
class ACSetMatcher(ExactMatcher):
	def __init__(self, patterns):
		# do initialization here 
		# don't use super call since that uses preprocessPattern which we don't want to do here
		self.preprocessPatterns(patterns)
		
	def findFailureLink(self,node,x):
		# implement algorithm of section 3.4.5
		# consult KMPMatcher to see how it uses the automata representation 
		# in that case
		pass

	def preprocessPatterns(self, patterns):
		# use findFailureLink here to construct complete keyword tree for
		# pattern set
		pass

	def matchTarget(self, target):
		# For testing purposes, note that this method as defined by ExactMatcher
		# implements the AC Search Algorithm of Section 3.4.4
		#
		# However, for the project you need to implement full AC Search Algorithm of Section 3.4.6
		pass