from proj2Lib.KMPMatcher import KMPMatcher as KMPMatcher

# UPDATE THIS CLASS TO IMPLEMENT SET MATCHING WITH THE KMP ALGORITHM
class NaiveSetMatcher(object):
	def __init__(self, patterns):
		pass

	def matchTarget(self, target):
		pass
		
	def getTotalTime(self):
		sum = 0
		for matcher in self.matchers:
			sum += matcher.getTotalTime()
		return sum
