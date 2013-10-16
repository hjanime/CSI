import numpy as np
#import getSummit as gs
import wig as WIG
import argparse
from operator import itemgetter
from sortedcollection import SortedCollection
import bisect
from dpgmm import dpgmm
import pylab as pl
from multiprocessing import Process, Queue, freeze_support

class Peak:
    '''
    A class to represent peaks with various functions such as dissect peaks etc.
    '''
    def __init__(self, ptuple):
        '''
        The coordinates will be 1-based.
        '''
        self.chrom, self.start, self.end, self.name, self.score, self.strand = ptuple
        self.dpgmmModel = None
        self.dpgmmDraw = None
        self.subPeaks = []
        self.totalTags = None
        self.wig = None
        self.color_iter = ['r','g','b','y','c','m']
        self.frac = 0.2
        self.selected = None
        self.dividePoints = None
        self.sbools = None
        self.parent = None
        self.density = None

    def setWig(self, chromWig):
        '''
        Set the original wig data within the peak region.
        Argument: chromWig, the wig data for the whole chromosome.
        Result: the self.wig will be set with the same columns as the chromWig.
        '''
        ws = bisect.bisect_left(chromWig[:,0], self.start)
        we = bisect.bisect_right(chromWig[:,0], self.end)
        self.wig = chromWig[ws:we,:]
        self.end = int(self.wig[-1,0]) + 1
        self.start = int(self.wig[0,0])
        self.totalTags = self.wig[:,1].sum()

    def __add_data_to_model(self):
        '''
        Helper method to populate the DPGMM model with the wig data.
        '''
        for i in range(self.wig.shape[0]):
            for j in range(int(round(self.wig[i,1]))):
                self.dpgmmModel.add([ self.wig[i,0] ])

    def dissect(self, n_components=6):
        '''
        Apply DPGMM to the peak and find the components.
        Argument: n_components (optional), the maximum number of components. Defaults to 5, which is sufficient in most cases.
        '''
        #Generate the DPGMM model
        self.dpgmmModel = dpgmm.DPGMM(1, n_components)
        self.__add_data_to_model()
        self.dpgmmModel.setPrior()
        self.dpgmmModel.solve()
        self.dpgmmDraw = np.array(self.dpgmmModel.intMixture())
        #Divide the peak based on the model.
        self._generateSelectedRegions()

        #Populate the subpeaks by the dividing points
        self.populateSubPeaks()

    def _generateSelectedBools(self):
        self.sbools = self.dpgmmDraw[0] > self.dpgmmDraw[0].max() * self.frac

    def _generateSelectedComps(self):
        '''
        Get the selected components.
        '''

        if self.sbools == None:
            self._generateSelectedBools()

        sweights = self.dpgmmDraw[0][self.sbools]

        scaleFactor = 1.0 / sweights.sum()
        '''
        print scaleFactor
        print sbools
        print self.dpgmmDraw[0][sbools]
        print self.dpgmmDraw[1][sbools]
        '''
        selected = zip(scaleFactor * self.dpgmmDraw[0][self.sbools], self.dpgmmDraw[1][self.sbools])
        selected.sort(key = lambda k: k[1].getLoc() )

        self.selected = selected

    def _generateSelectedRegions(self):
        '''
        Find the sub-regions that are selected and belong to a component.
        The self.subRegions should be in the format (start, end, component_idx)
        '''
        if self.sbools == None:
            self._generateSelectedBools()
        if self.selected == None:
            self._generateSelectedComps()
        axis_x = range(int(self.wig[0,0]), int(self.wig[-1,0]) + 1)
        stick = []
        for x in axis_x:
            stick.append( np.argmax(self.dpgmmModel.stickProb([x])) )
        subRegions = []
        '''
        lastStick = stick[0]
        lastPos = 0
        for i,l in enumerate(stick):
            if l != lastStick:
                if self.sbools[lastStick]:
                    subRegions.append((lastPos, i-1, lastStick))
                lastPos = i
                lastStick = l
        if self.sbools[lastStick]:
            subRegions.append((lastPos, i-1, lastStick))
        '''

        stick = np.array(stick)
        maxs, mins = self._getLocalMaxAndMin(self.selected)
        lastMin = 0
        mins = list(mins)
        mins.append(len(axis_x))
        for i in list(mins):
            startPos = lastMin
            endPos = i
            for j in range(lastMin, i):
                if self.sbools[stick[j]]:
                    startPos = j
                    break
            j = i - 1
            while not self.sbools[stick[j]] and j >= startPos:
                j -= 1
            endPos = j + 1
            if endPos - startPos > 1:
                #The commented code is used to find the dominant component
                #tempCounts = np.bincount(stick[startPos:endPos])
                #maxStick = np.argmax(tempCounts)

                subRegions.append((startPos, endPos, list(set(stick[startPos:endPos]))))
            lastMin = i



        self.subRegions = subRegions


    def setParent(self, peak):
        self.parent = peak

    def populateSubPeaks(self):
        '''
        Calculate the subpeaks given the mixture model.
        '''
        '''
        if not self.selected:
            self._getSelectedComps()
        selected = self.selected
        maxs, mins = self._getLocalMaxAndMin(selected)
        '''
        self._generateDensity()
        if len(self.subRegions) > 1:
            for r in self.subRegions:
                self.subPeaks.append(Peak((self.chrom, self.start + r[0], self.start + r[1], self.name+'_1', self.score, self.strand  )))
                self.subPeaks[-1].setParent(self)
                ws = bisect.bisect_left(self.wig[:,0], self.start+r[0])
                we = bisect.bisect_right(self.wig[:,0],self.start + r[1])
                self.subPeaks[-1].setWig(self.wig[ws:we,:])
            #print self.asBedStringRe()
            return True
        else:
            return False


    def _getProb(self, selected,x):
        '''
        Return the probability density of a point.
        '''
        y = 0
        for s in selected:
            #print s
            y += s[0] * s[1].prob(x)

        return y

    def _getComponentProbs(self):
        '''
        Return a array of the components' probability on all points within the region.
        '''
        if self.selected == None:
            self._generateSelectedComps()
        axis_x = range(int(self.wig[0,0]), int(self.wig[-1,0]) + 1)
        y = []
        for s in self.selected:
            tempY = []
            for x in axis_x:
                tempY.append( s[0] * s[1].prob(x) )
            y.append(tempY)
        y= np.array(y)
        return np.array(axis_x), y

    def _generateDensity(self):
        x, ys = self._getComponentProbs()
        totalReads = self.wig[:,1].sum()
        cy = ys.sum(0)
        self.density = cy*totalReads

    def getDensity(self):
        if self.density != None:
            return self.density
        elif self.density == None and self.parent == None:
            self._generateDensity()
            return self.density
        elif self.density == None:
            assert self.parent.density != None
            return self.parent.density[self.start - self.parent.start : self.end - self.parent.start]

    def getScoreByDensity(self, apeak):
        '''
        Return the pairing score of this peak and another one by using density.
        '''
        assert self.strand != apeak.strand
        if self.strand == '-':
            return apeak.getScoreByDensity(self)
        density1 = self.getDensity()
        density2 = apeak.getDensity()
        #bestScore = sys.maxint
        bestScore = 0
        bestShift = 0
        bestHeight = 0
        bestAve = 0
        bestRpos = 0
        shift = 0
        while self.start + shift < apeak.end:
            start1 = self.start + shift
            currStart = max(start1, apeak.start)
            currEnd = min(self.end + shift, apeak.end)
            if currStart >= currEnd:
                shift += 1
                continue
            reg1 = density1[currStart - start1: currEnd - start1 ]
            reg2 = density2[currStart - apeak.start: currEnd - apeak.start ]
            score = 0
            tempHeight = 0
            total = 0
            tempRpos = 0
            for idx, (i, j) in enumerate(zip(reg1, reg2)):
                score += i*j
                if i * j > tempHeight:
                    tempHeight = i * j
                    tempRpos = currStart + idx
                total += i * j
            tempAve = total / (currEnd - currStart)

            #score += currStart - min(start1, apeak.start) + max(self.end+shift, apeak.end) - currEnd
            if score > bestScore:
                bestShift = shift
                bestScore = score
                bestHeight = tempHeight
                bestAve = tempAve
                bestRpos = tempRpos
            shift += 1
        return bestScore, bestShift, bestRpos, bestAve, bestHeight

    def getCombinedDensity(self, apeak, shift):
        '''
        Get the combined density profiles of the current peak and its mate.
        '''
        if self.strand == '-':
            return apeak.getCombinedDensity(self)
        cstart = max(self.start + shift, apeak.start)
        cend = min(self.end + shift, apeak.end)
        fregion = self.getDensity()[cstart - self.start - shift : cend - self.start - shift]
        rregion = apeak.getDensity()[cstart - apeak.start : cend - apeak.start]
        combined = np.zeros(cend-cstart)
        try:
            combined = fregion + rregion
        except ValueError:
            print "shift: ", shift
            print self.density.shape, ' ', apeak.density.shape
            print self.asBedString().strip()
            print apeak.asBedString()
        start = cstart - shift/2
        end = cend - shift/2
        return (start, end, combined)


    def plotCompSelected(self,savefig=True):
        '''
        Plot the seleceted components.
        '''
        x, ys = self._getComponentProbs()
        totalReads = self.wig[:,1].sum()
        pl.scatter(self.wig[:,0], self.wig[:,1], color='0.5')
        for y,c in zip(ys, self.color_iter):
            pl.scatter(x, y*totalReads, color=c)

        cy = ys.sum(0)
        for dp,c in zip(self.subRegions, self.color_iter):
            pl.scatter(x[dp[0]:dp[1]+1], cy[dp[0]:dp[1]+1]*totalReads, color=c)

        if savefig:
            #pl.scatter(x[lastDp:], cy[lastDp:]*totalReads*5, color='k')
            outDir = 'test_png_smoothed_nosub_lowlim100/' + self.chrom + '/strand_' + self.strand + "/"
            import os
            if not os.path.exists(outDir):
                os.makedirs(outDir)

            pl.savefig(outDir + '_'.join([self.name.replace('/','_'), self.chrom, str(self.start), str(self.end)])+ '.png')
            pl.clf()

    def plotShiftedAll(self, rpeak, shift, savefig = True):
        assert self.strand != rpeak.strand
        if self.strand == '-':
            rpeak.plotShiftedAll(self, shift)
        else:
            combinedDensity = self.getCombinedDensity(rpeak, shift)
            pl.scatter(range(self.start + shift/2, self.end + shift/2), self.getDensity(), color="r")
            pl.scatter(range(rpeak.start - shift/2, rpeak.end - shift/2), rpeak.getDensity(), color="b")
            pl.scatter(range(combinedDensity[0], combinedDensity[1]), combinedDensity[2], color="g")
            pl.scatter(self.wig[:,0]+shift/2, self.wig[:,1],color='k')
            pl.scatter(rpeak.wig[:,0]-shift/2, rpeak.wig[:,1], color='0.5')
            if savefig:
                outDir = 'paired_png_smoothed_nosub_lowlim100/' + self.chrom + '/'
                import os
                if not os.path.exists(outDir):
                    os.makedirs(outDir)
                pl.savefig(outDir + '_'.join([self.name.replace('/','_'), rpeak.name.replace('/','_'), self.chrom, str(shift), str(self.start)]) + '.png')
                pl.clf()



    def _getLocalMaxAndMin(self, selected):
        '''
        Return two lists, one for local maximums and one for local minimums.
        The local minimums should be located between local maximums, as dividing points.
        '''
        x = range(int(self.wig[0,0]), int(self.wig[-1,0]) + 1)
        y = []
        for i in x:
            y.append(self._getProb(selected, i))

        from scipy.signal import argrelextrema
        y = np.array(y)
        #print y
        maxs = argrelextrema(y, np.greater)[0]
        mins = argrelextrema(y, np.less)[0]
        #print mins
        #print maxs

        return maxs, mins

    def asBedStringRe(self):
        '''
        Return a string of the peaks in BED format.
        For each peak with subpeaks, a hierarchical method will be used.
        A peak will be followed by all its subpeaks ordered by location.
        '''
        strs = []
        strs.append(self.asBedString())
        for sp in self.subPeaks:
            strs.append(sp.asBedString())
        return ''.join(strs)


    def asBedString(self):
        '''
        Return a string of the peak in BED format.
        Only including the current peak.
        '''
        return '%s\t%d\t%d\t%s\t%f\t%s\n'%(self.chrom, self.start-1, self.end, self.name, self.score, self.strand)

    def asList(self):
        '''
        Return a list of the data.
        '''
        return [self.chrom, self.start-1, self.end, self.name, self.score, self.strand]











def loadPeaks( filename ):
    f = open( filename )
    fpeaks = {}
    rpeaks = {}
    for r in f:
        if not ( r.startswith('#') or r.startswith("track") ):
            tokens = r.strip().split()
            chrom = tokens[0]
            start = int( tokens[1] ) + 1 #Adjust to 1-based coordinates.
            end = int( tokens[2] )
            name = ''
            score = 0
            strand = '.'
            if len( tokens ) > 3:
                name = tokens[3]
            if len( tokens ) > 4:
                score = float( tokens[4] )
            if len( tokens ) > 5:
                strand = tokens[5]
            newPeak = Peak((tokens[0], start, end, name, score, strand))
            #newPeak.dissect()
            if strand == '-':
                if chrom not in rpeaks:
                    rpeaks[ chrom ] = []
                rpeaks[chrom].append( newPeak )
            else:
                if chrom not in fpeaks:
                    fpeaks[ chrom ] = []
                fpeaks[chrom].append( newPeak )
    f.close()

    for chrom in fpeaks:
        fpeaks[chrom].sort(key = lambda k:( k.chrom, k.start, k.end ))
        if chrom in rpeaks:
            rpeaks[chrom].sort(key = lambda k:( k.chrom, k.start, k.end ))
    return fpeaks, rpeaks

def getScore( fp, rp, fstart, rstart ):
    '''
    fp and rp are two expanded arrays of the tag counts.
    They are expanded as end-start+1 as in getSummit expandWig.

    Returns the score and the shift of the two peaks
    '''
    fend = fstart + len( fp ) - 1
    rend = rstart + len( rp ) - 1
    fmean = sum( fp ) / len( fp )
    rmean = sum( rp ) / len( rp )
    newfp = np.array( [ ( t / fmean) ** 2 for t in fp ], dtype="float64")
    newrp = np.array( [ ( t / rmean) ** 2 for t in rp ], dtype="float64")
    maxScore = 0
    shift = 0
    i = 0
    maxRpos = 0
    maxheight = 0
    average = 0
    #Gradually move forward towards the 3' end
    while fstart + i < rend:
        overlap = 0
        tempStart = max( fstart + i, rstart )
        tempEnd = min( fend + i, rend )
        tempRpos = 0
        tempHeight = 0
        totalRaw = 0
        for j in range( tempStart, tempEnd + 1 ):
            temp =  newfp[ j - ( fstart + i ) ] * newrp[ j - rstart ]
            tempRawHeight = fp[ j - ( fstart + i ) ] + rp[ j - rstart ]
            if tempRawHeight >  tempHeight:
                tempHeight = tempRawHeight
                tempRpos = j
            overlap += temp
            totalRaw += fp[ j - ( fstart + i ) ] + rp[ j - rstart ]
        tempScore =  overlap #/ ( abs (sum( rp[ tempStart - rstart : tempEnd - rstart + 1] ) - sum( fp[ tempStart - (fstart+i) : tempEnd - (fstart+i) +1] )) + 10 )
        if tempScore > maxScore:
            maxScore = tempScore
            shift = i
            maxheight = tempHeight
            maxRpos = tempRpos
            average = totalRaw / ( tempEnd - tempStart + 1 )
        i += 1

    return maxScore, shift, maxRpos, average, maxheight



def pair( fpeaks, rpeaks, ulimit, dlimit):
    '''
    Assuming that the peaks on one strand is mutually exclusive.
    They do not overlap with each other.
    In this case, the ordering of the starts of the peaks and the
    ends of the peaks are the same. And that when the starts are
    sorted, the ends are also sorted.
    '''
    #print "fwig: ", fwig
    #print "rwig: ", rwig
    #offset = 5
    #expandCol = 1
    #out1 = open(prefix + "_singletons_shape_c2c.bed",'w')
    #out2 = open(prefix + "_pairs_shape_c2c.gff", "w") #c2c means center-to-center
    #out3 = open(prefix + "_pairs_shape_true_e2e_detail.gff",'w') #e2e means end-to-end.
    #out4 = open(prefix + "_pairs_shape_true_e2e.gff",'w')
    #print fpeaks.keys()
    #for chrom in fpeaks:
    #    if chrom not in rpeaks:
    #        continue
    #    print chrom
    fp = fpeaks
    rp = rpeaks
    pairF = []  #Store the pairing information, if unpaired, it will be negative.
    pairR = []
    #fw = fwig[chrom]
    #expandedFw = WIG.expandWig( fw, offset, expandCol, smooth=False )
    #rw = rwig[chrom]
    #expandedRw = WIG.expandWig( rw, offset, expandCol, smooth=False )
    rstarts = []
    rends = []
    fprefer = []
    rprefer = []
    unpairedF = []
    for f in fp:
        fprefer.append( SortedCollection( key=itemgetter(1) ) )
        pairF.append( ( -1, 0, 0) )  #( index of the mate, score, distance )
    for r in rp:
        rprefer.append( SortedCollection( key=itemgetter(1) ) )
        pairR.append( (-1, 0, 0) )
        rstarts.append( r.start )
        rends.append( r.end )
    for i in range( len( fp ) ):
        currfp = fp[ i ]
        start = currfp.start
        end = currfp.end
        es = start - ulimit
        ee = end + dlimit
        #currFw = expandedFw[ max( 0, start - fw[0,0] ) + offset : max( 0, end - fw[0,0] ) + offset + 1]
        si = bisect.bisect_left( rends, es )
        ei = bisect.bisect_right( rstarts, ee )
        #ftagCounts,_,_ = gs.getTagCount( fwig, chrom, start, end )

        #print ei - si
        maxScore = 0
        bestDist = 0
        bestIdx = 0
        bestRpos = 0
        bestAve = 0
        bestHeight = 0
        for idx in range( si, ei ):
            currrp = rp[ idx ]
            #rstart = currrp.start
            #rend = currrp.end
            #currRw = expandedRw[ max( 0, rstart - rw[0, 0] ) + offset : max( 0, rend - rw[0,0] ) + offset + 1 ]
            #rtagCoungs,_,_ = gs.getTagCount( rwig, chrom, currrp[1], currrp[2] )

            tempScore, tempDist, tempRpos, tempAve, tempHeight = currfp.getScoreByDensity( currrp )

            fprefer[ i ].insert( (idx, tempScore, tempDist, tempRpos, tempAve, tempHeight) )
            rprefer[ idx ].insert( (i, tempScore, tempDist, tempRpos, tempAve, tempHeight) )
            if tempScore > maxScore:
                maxScore = tempScore
                bestDist = tempDist
                bestIdx = idx
                bestRpos = tempRpos
                bestAve = tempAve
                bestHeight = tempHeight
        if maxScore > pairR[ bestIdx ][1]:
            pairF[i] = ( bestIdx, maxScore, bestDist , bestRpos, bestAve, bestHeight)
            if pairR[ bestIdx ][0] >= 0:
                pairF[ pairR[ bestIdx ][ 0 ] ] = (-1, 0, 0)
                unpairedF.append( pairR[ bestIdx ][0] )
            pairR[bestIdx] = ( i, maxScore, bestDist, bestRpos, bestAve, bestHeight )
        else:
            unpairedF.append( i )
        try:
            fprefer[ i ].remove( ( bestIdx, maxScore, bestDist, bestRpos, bestAve, bestHeight ) )
        except ValueError:
            #print "Value error: ", bestIdx, ' ',maxScore, ' ',bestDist,' ', si,' ', ei
            pass
    singletons = []
    pairs = []
    pairs_extended = []
    combinedDensity = []
    while len(unpairedF) > 0:
        for u in unpairedF:
            if len( fprefer[u] ) > 0:
                ridx = fprefer[u][-1][0]
                if pairR[ ridx ][1] < fprefer[u][-1][1]:
                    if pairR[ ridx ][0] > 0:
                        pairF[ pairR[ ridx ][0] ] = (-1, 0, 0)
                        unpairedF.append( pairR[ ridx ][0] )
                    pairR[ ridx ] = ( u, fprefer[u][-1][1], fprefer[u][-1][2], fprefer[u][-1][3], fprefer[u][-1][4], fprefer[u][-1][5] )
                    pairF[ u ] = ( ridx, fprefer[u][-1][1], fprefer[u][-1][2], fprefer[u][-1][3], fprefer[u][-1][4], fprefer[u][-1][5] )
                fprefer[u].remove( (ridx, fprefer[u][-1][1], fprefer[u][-1][2], fprefer[u][-1][3], fprefer[u][-1][4], fprefer[u][-1][5] ) )
            else:
                unpairedF.remove( u )

    for i,f in enumerate(pairF):
        #fp[i][1] -= 1
        if f[0] == -1:
            singletons.append( fp[i] )
        else:
            #rp[f[0]][1] -= 1
            #pairs.append( fp[i] )
            #pairs.append( rp[f[0]] )
            #pairStart = (2*f[3] - f[2])/2
            #pairEnd = pairStart + 1
            pairStart = f[3] - f[2] + 1
            pairEnd = f[3]
            pairs.append( [fp[i].chrom,f[5],'.',pairStart, pairEnd,f[4],'.','.','cw_distance='+str(f[2]) ] )
            pairs_extended.append( [fp[i].chrom,'.','.',fp[i].start+1, rp[f[0]].end,f[5],'.','.','bestAve=%f,maxScore=%f,bestDist=%d'%(f[4],f[1],f[2]) ] + fp[i].asList() + rp[f[0]].asList() )
            combinedDensity.append(fp[i].getCombinedDensity(rp[f[0]], f[2]))
            fp[i].plotShiftedAll(rp[f[0]], f[2])

    for i,f in enumerate(pairR):
        if f[0] == -1:
            singletons.append( rp[i] )

    singletons.sort(key=lambda k:( k.chrom, k.start, k.end ))
    pairs.sort(key = lambda k:( k[0], k[1], k[2]))
    print "singletons: ", len(singletons)
    print "pairs: ", len(pairs)

    return singletons, pairs, pairs_extended, combinedDensity


def pair_region(fpeaks, rpeaks, extend=200):
    '''
    1. Dissect both forward and reverse peaks.
    2. For each subpeak of the peaks on the forward strand, find the matching scores of each subregion on the
       reverse strand and rank them by matching score from high to low.
    3. Use dynamic programming to find the matches that gives the best overall score for the region with enforcement
       of the order.
    '''
    for fp in fpeaks:
        subPeaks = []
        if fp.subPeaks != None and len(fp.subPeaks) > 1:
            for sp in fp.subPeaks:
                subPeaks.append(sp)
        else:
            subPeaks.append(fp)


def processChromPeaks(taskQ, outQ, fwig, rwig, ulimit, dlimit, processID):
    for chrom, fpeaks, rpeaks in iter(taskQ.get, "STOP"):
        print "Process ", processID, ' is processing peaks in ', chrom
        allfpeaks = []
        allrpeaks = []
        for p in fpeaks:
            p.setWig(fwig[chrom])
            p.dissect()
            p.plotCompSelected()
            #outQ.put(p.asBedStringRe())
            '''
            if p.subPeaks != None and len(p.subPeaks) > 1:
                allfpeaks += p.subPeaks
            else:
                allfpeaks.append(p)
            '''
            allfpeaks.append(p)
        for p in rpeaks:
            p.setWig(rwig[chrom])
            p.dissect()
            p.plotCompSelected()
            '''
            if p.subPeaks != None and len(p.subPeaks) > 1:
                allrpeaks += p.subPeaks
            else:
                allrpeaks.append(p)
            '''
            allrpeaks.append(p)
        print "Process ", processID, " is pairing"
        output = pair(allfpeaks, allrpeaks, ulimit, dlimit)
        outQ.put((chrom, output))
        print "Process ", processID, " finished."

def writePeaks(peaks, rpeaks, filename,fwig, rwig, ulimit, dlimit, prefix):
    '''
    Write the peaks in bed format.
    '''
    freeze_support()
    outQ = Queue()
    NUM_PROCESSES = 3
    taskQ = Queue()
    count = 0
    for chrom in peaks:
        if chrom in rpeaks:
            taskQ.put((chrom, peaks[chrom], rpeaks[chrom]))
            count += len(peaks[chrom]) + len(rpeaks[chrom])

    processes = []
    for i in range( NUM_PROCESSES ):
        processes.append(Process( target=processChromPeaks, args=(taskQ, outQ, fwig,rwig, ulimit, dlimit, i)))
        processes[-1].start()

    out1 = open(prefix + "_singletons_shape_c2c.bed",'w')
    out2 = open(prefix + "_pairs_shape_c2c.gff", "w") #c2c means center-to-center
    out3 = open(prefix + "_pairs_shape_true_e2e_detail.gff",'w') #e2e means end-to-end.
    out4 = open(prefix + "_pairs_shape_true_e2e.gff",'w')
    out5 = open(prefix + "_combined_density.bedgraph", "w")
    #f = open(filename, 'w')
    '''
    for chrom in peaks:
        for p in peaks[chrom]:
            p.setWig(wig[chrom])
            p.dissect()
            p.plotCompSelected()
            f.write(p.asBedStringRe())
    f.close()
    '''

    for i in range(count):
        chrom,(singletons, pairs, pairs_extended,combinedDensity) =  outQ.get()

        for s in singletons:
            out1.write(s.asBedString())
        for p in pairs:
            out2.write('\t'.join([str(i) for i in p]))
            out2.write('\n')
        out3.write('#')
        out3.write('\t'.join(['chrom','source','feature','start','end','bestHeight','frame','aattribute','peak1','','','','','', 'peak2']))
        out3.write('\n')
        for p in pairs_extended:
            out3.write('\t'.join([str(i) for i in p]))
            out3.write('\n')
            out4.write('\t'.join([str(i) for i in p[0:8]]))
            out4.write('\n')
        for s, e, values in combinedDensity:
            for i in range(s, e):
                out5.write('\t'.join([chrom, str(i), str(i+1), str(values[i-s])]))
                out5.write('\n')
    out1.close()
    out2.close()
    out3.close()
    out4.close()
    out5.close()


    #f.close()
    for i in range( NUM_PROCESSES ):
        taskQ.put("STOP")
    for i in range( NUM_PROCESSES ):
        processes[i].join()


def main():
    parser = argparse.ArgumentParser(description="Pair the peaks on different strand.")
    parser.add_argument( "peak", help="The bed file contains all peaks." )
    parser.add_argument( 'forwardWig', help="The forward wig file." )
    parser.add_argument( 'reverseWig', help="The reverse wig file." )
    parser.add_argument( 'output', help="The output prefix. Will output two files, one containing all the singletons and the other with the pairs. The pair file will be in GFF format and singletons in BED format. Output files will be prefix_singletons.bed and prefix_pairs.gff" )
    parser.add_argument( '-d', '--downstream', type = int, default = 100, help = "Within how far downstream to look for a mate." )
    parser.add_argument( '-u', '--upstream', type = int, default = 0, help = "Within how far upstream to look for a mate." )

    args = parser.parse_args()

    fpeaks, rpeaks = loadPeaks( args.peak )
    fwig = WIG.loadWig( args.forwardWig, smooth=False )
    rwig = WIG.loadWig( args.reverseWig, smooth=False, strand='-' )
    writePeaks( fpeaks, rpeaks, 'test_peaks_smoothed_lowlim100.bed', fwig, rwig, args.upstream, args.downstream, args.output)

    #pair( fpeaks, rpeaks, fwig, rwig, args.upstream, args.downstream, args.output)

if __name__ == "__main__":
    main()
