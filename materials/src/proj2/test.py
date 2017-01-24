from proj2Lib.NaiveMatcher import NaiveMatcher as NaiveMatcher
from proj2Lib.KMPMatcher import KMPMatcher as KMPMatcher

pattern = 'abbabbabb'

target = 'abbabbabba'

for i in xrange(10):
	target += target

m1 = NaiveMatcher(pattern)
print m1

m1_ptime = m1.getTotalTime()

m1.matchTarget(target)
m1_mtime = m1.getTotalTime() - m1_ptime


m2 = KMPMatcher(pattern)
print m2

m2_ptime = m2.getTotalTime()

m2.matchTarget(target)
m2_mtime = m2.getTotalTime() - m2_ptime

print 'm1: preproc: %d, match: %d' % (m1_ptime, m1_mtime)
print 'm2: preproc: %d, match: %d' % (m2_ptime, m2_mtime)


