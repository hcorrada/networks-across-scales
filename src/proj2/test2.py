from proj2Lib.NaiveSetMatcher import NaiveSetMatcher
from proj2Lib.KMPSetMatcher import KMPSetMatcher
from proj2Lib.ACSetMatcher import ACSetMatcher

patterns = ['potato', 'tattoo', 'theater', 'other']
target = 'xxpotattooxx'

nsm = NaiveSetMatcher(patterns)
nsm_ptime = nsm.getTotalTime()

nsm.matchTarget(target)
nsm_mtime = nsm.getTotalTime() - nsm_ptime

print 'nsm: ptime %d, mtime %d' % (nsm_ptime, nsm_mtime)

ksm = KMPSetMatcher(patterns)
ksm_ptime = ksm.getTotalTime()

ksm.matchTarget(target)
ksm_mtime = ksm.getTotalTime() - ksm_ptime

print 'ksm: ptime %d, mtime %d' % (ksm_ptime, ksm_mtime)

acm = ACSetMatcher(patterns)
acm_ptime = acm.getTotalTime()

print acm
acm.matchTarget(target)
acm_mtime = acm.getTotalTime() - acm_ptime
print 'acm: ptime %d, mtime %d' % (acm_ptime, acm_mtime)

tes=ACSetMatcher(['potato', 'tatter', 'at'])
print tes

tes.matchTarget('potatter')

tes=ACSetMatcher(['acatt', 'ca'])
print tes
tes.matchTarget('acatg')

