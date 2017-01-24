from proj2Lib.NaiveMatcher import NaiveMatcher as NaiveMatcher

# a Naive Pattern Set exact matcher
# uses the naive exact matching algorithm for each pattern and the target
# use this as basis for your KMPSetMatcher implementation
class NaiveSetMatcher(object):
	def __init__(self, patterns):
		self.matchers = [NaiveMatcher(patterns[i], reportFormat = 'Found match of pattern %d at position {0}' % i) for i in xrange(len(patterns))]

	def matchTarget(self, target):
		for matcher in self.matchers:
			matcher.matchTarget(target)

	def getTotalTime(self):
		sum = 0
		for matcher in self.matchers:
			sum += matcher.getTotalTime()
		return sum
