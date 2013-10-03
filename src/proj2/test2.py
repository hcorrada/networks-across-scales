from proj2Lib.NaiveSetMatcher import NaiveSetMatcher as NaiveSetMatcher

patterns = ['potato', 'tattoo', 'theater', 'other']
target = 'xxpotattooxx'

nsm = NaiveSetMatcher(patterns)
nsm_ptime = nsm.getTotalTime()

nsm.matchTarget(target)
nsm_mtime = nsm.getTotalTime() - nsm_ptime

print 'nsm: ptime %d, mtime %d' % (nsm_ptime, nsm_mtime)