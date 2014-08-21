import numpy as np

def preprocMat(scoringMatrix):
    keys = scoringMatrix.keys()
    alphabet = set(zip(*keys)[0])
    mat = np.zeros([len(alphabet), len(alphabet)])
    alphabetInd = dict(zip(alphabet, range(len(alphabet))))
    for ((x,y), val) in scoringMatrix.iteritems():
        xind = alphabetInd[x]
        yind = alphabetInd[y]
        mat[xind,yind] = val
        mat[yind,xind] = val
    return mat, alphabetInd

class AffineScoring(object):
    def __init__(self, scoringMatrix, gapStart, gapExtend):
        self.scoringMatrix = scoringMatrix
        self.mat, self.alphabetInd = preprocMat(scoringMatrix)
        self.gapStart = gapStart
        self.gapExtend = gapExtend

    def match(self, tup):
        if tup not in self.scoringMatrix:
            tup = tup[::-1]
        return self.scoringMatrix[tup]

    def start(self):
        return self.gapStart

    def extend(self):
        return self.gapExtend

    def indexify(self,x,y):
        xind = [self.alphabetInd[xc] for xc in x]
        yind = [self.alphabetInd[yc] for yc in y]
        return xind, yind
    
    def matchMat(self,xcind,yind):
        return self.mat[xcind,yind].reshape([len(yind),1])

class LinearScoring(AffineScoring):
    def __init__(self, scoringMatrix, gapPenalty):
        super(LinearScoring, self).__init__(scoringMatrix, 0, gapPenalty)

class ArrowKind:
    none = 0
    m = 1
    ix = 2
    iy = 3
    
class Alignment(object):
    def __init__(self, x, y, trace, score, nIdentities, nGaps):
        self.trace = trace
        self.score = score
        self.nIdentities = nIdentities
        self.nGaps = nGaps
        self.x = x
        self.y = y
        self.nx = len(x)
        self.ny = len(y)
        
    def __len__(self):
        return len(self.trace)

    @property
    def alnLen(self):
        return len(self.trace)
    
    @property
    def propId(self):
        return float(self.nIdentities) / self.alnLen

    @property
    def propGap(self):
        return float(self.nGaps)/ self.alnLen 

    def getHeader(self):
        return 'Score = {0.score}, Identities = {0.nIdentities}/{0.alnLen} ({0.propId:.0%}), Gaps = {0.nGaps}/{0.alnLen} ({0.propGap:.0%})\n'.format(self)

    def getAlnStrings(self, start, stop):
        alnx = ''
        alny = ''
        alnm = ''

        for k in xrange(start, stop):
            arrow = self.trace[k]
            transKind = arrow[0]
            i,j = arrow[1]

            xchar = self.x[i-1]
            ychar = self.y[j-1]
            
            if transKind == ArrowKind.ix:
                alnx += '-'
                alny += ychar
                alnm += ' '
                
            elif transKind == ArrowKind.iy:
                alnx += xchar
                alny += '-'
                alnm += ' '
            else:
                alnx += xchar
                alny += ychar
                alnm += xchar if xchar == ychar else ' '
        return alnx, alny, alnm

    def getPartialAln(self, start, stop):
        alnx, alny, alnm = self.getAlnStrings(start, stop)
        startLoc = self.trace[start][1]
        endLoc = self.trace[stop-1][1]

        out = 'Query\t{0}\t{1}\t{2}\n'.format(startLoc[0], alnx, endLoc[0])
        out += '\t\t{0}\n'.format(alnm)
        out += 'Text\t{0}\t{1}\t{2}\n'.format(startLoc[1], alny, endLoc[1])
        return out
    
    def getAln(self, width=80):
        out = self.getHeader()
        out += '\n'

        start = 0
        while start <= len(self):
            stop = start + width
            out += self.getPartialAln(start, min(stop, len(self))) + '\n'
            start = stop
        return out

    def __repr__(self):
        return self.getAln(width=80)

class AlnMatrix(object):
    def __init__(self, nx, ny):
        self._val = np.zeros([nx+1,ny+1,4], dtype=np.float)
        self._arrow = np.zeros([nx+1,ny+1,4], dtype=np.int)

    @property
    def val(self):
        return self._val

    @property
    def arrow(self):
        return self._arrow
    
    @property
    def shape(self):
        return self._val.shape
    
    def reprCell(self, i, j, arrowKind):
        return '%.0f:%d' % (self.val[i,j,arrowKind], self.arrow[i,j, arrowKind])

    def rowStrings(self, arrowKind):
        nx, ny, _ignore = self.shape
        strs = ['' for j in xrange(ny)]
        for i in xrange(nx):
            for j in xrange(ny):
                strs[j] += ('\t' + self.reprCell(i, j, arrowKind))
        return strs

class Aligner(object):
    def __init__(self, x, y, score, alnKind):
        self.x = x
        self.y = y
        self.nx = len(x)
        self.ny = len(y)
        self.alnKind = alnKind
        self.score = score
        self._m = AlnMatrix(self.nx, self.ny)
        self.xind, self.yind = score.indexify(x,y)
        
    @property
    def val(self):
        return self._m.val

    @property
    def arrow(self):
        return self._m.arrow
        
    def reprMat(self, arrowKind):
        rowStrings = self._m.rowStrings(arrowKind)
        tmp = '*' + self.y
        for i in xrange(len(rowStrings)):
            rowStrings[i] = tmp[i] + rowStrings[i]
        tmp = '\t*\t' + '\t'.join(list(self.x))
        return tmp + '\n' + '\n'.join(rowStrings) + '\n'
    
    def __repr__(self):
        out =  'm matrix:' + '\n'
        out += self.reprMat(ArrowKind.m) + '\n'

        out +=  '\nix matrix:' + '\n'
        out += self.reprMat(ArrowKind.ix) + '\n'

        out += '\niy matrix:' + '\n'
        out += self.reprMat(ArrowKind.iy)

        return out

    def init(self):
        self.val[0,0,ArrowKind.ix] = -np.inf
        self.val[0,0,ArrowKind.iy] = -np.inf

        nx = self.nx
        ny = self.ny
        score = self.score 
        alnKind = self.alnKind
        
        # fill-in first rows
        self.val[1:,0,ArrowKind.m] = -np.inf
        self.val[1:,0,ArrowKind.ix] = -np.inf
        self.val[1:,0,ArrowKind.iy] = score.start() + score.extend() * np.arange(1,nx+1) if alnKind == 'global' else 0
        self.arrow[1:,0,ArrowKind.iy] = np.repeat(ArrowKind.iy if alnKind == 'global' else ArrowKind.none, nx)
        
        # fill-in first columns
        self.val[0,1:,ArrowKind.m] = -np.inf
        self.val[0,1:,ArrowKind.iy] = -np.inf
        self.val[0,1:,ArrowKind.ix] = score.start() + score.extend() * np.arange(1,ny+1) if alnKind == 'global' else 0
        self.arrow[0,1:,ArrowKind.ix] = np.repeat(ArrowKind.ix if alnKind == 'global' else ArrowKind.none, ny)

    def fill(self):
        nx = self.nx
        ny = self.ny
        score = self.score
        
        localval =  0 if self.alnKind == 'local' else -np.inf
        scorearray = localval * np.ones([self.ny,4])

        iyarray = np.zeros([self.ny,4])
        iyarray[:,ArrowKind.none] = localval
        iyarray[:,ArrowKind.m] = score.start() + score.extend()
        iyarray[:,ArrowKind.ix] = score.start() + score.extend()
        iyarray[:,ArrowKind.iy] = score.extend()

        ixarray = np.zeros([1,4])
        ixarray[:,ArrowKind.none] = localval
        ixarray[:,ArrowKind.m] = score.start() + score.extend()
        ixarray[:,ArrowKind.ix] = score.extend()
        ixarray[:,ArrowKind.iy] = score.start() + score.extend()
        
        for i in xrange(1,nx+1):
            scorearray[:,1:] = score.matchMat(self.xind[i-1], self.yind)
            mchoices = scorearray + self.val[i-1,:-1,:]
            bestKind = np.argmax(mchoices, axis=1)
            self.arrow[i,1:,ArrowKind.m] = bestKind
            self.val[i,1:,ArrowKind.m] = mchoices[np.arange(ny), bestKind]

            iychoices = iyarray + self.val[i-1,1:,:]
            bestKind = np.argmax(iychoices, axis=1)
            self.arrow[i,1:,ArrowKind.iy] = bestKind
            self.val[i,1:,ArrowKind.iy] = iychoices[np.arange(ny), bestKind]

            for j in xrange(1,ny+1):
                ixchoices = ixarray + self.val[i,j-1,:]
                bestKind = np.argmax(ixchoices, axis=1)
                self.arrow[i,j,ArrowKind.ix] = bestKind
                self.val[i,j,ArrowKind.ix] = ixchoices[0, bestKind]
            
    def backtrace(self):
        nx = self.nx
        ny = self.ny
        
        if self.alnKind == 'global':
            bestKind = np.argmax(self.val[nx,ny,1:])
            bestArrow = (bestKind + 1, (nx,ny))
        else:
            i,j,bestKind = np.unravel_index(np.argmax(self.val[:,:,1:]), self.val[:,:,1:].shape)
            bestArrow = (bestKind + 1, (i,j))

        kind, (i,j) = bestArrow
        score = self.val[i,j,kind]
        nIdentities = 0
        nGaps = 0
        
        trace = []
        
        arrow = bestArrow        
        transKind = arrow[0]
        while (transKind != ArrowKind.none):
            # check if we're pointing to 0,0
            if i+j == 0:
                break

            # stop if this is a local alignment and we are pointing to a 0 position
            if self.alnKind == 'local' and self.val[i,j,transKind] == 0:
                break
                
            trace.append(arrow)
            nIdentities += 1 if transKind == ArrowKind.m and self.x[i-1] == self.y[j-1] else 0
            nGaps += 1 if transKind == ArrowKind.ix or transKind == ArrowKind.iy else 0

            # get the next arrow
            kind = self.arrow[i,j,transKind]
            if transKind in [ArrowKind.m, ArrowKind.iy]:
                i -= 1
            if transKind in [ArrowKind.m, ArrowKind.ix]:
                j -= 1

            arrow = (kind, (i,j))
            transKind = arrow[0]
        return Alignment(self.x, self.y, trace[::-1], score, nIdentities, nGaps)
    
def align(x, y, score, alnKind='global'):
    mat = Aligner(x, y, score, alnKind)
    mat.init()
    mat.fill()
    return mat.backtrace()
