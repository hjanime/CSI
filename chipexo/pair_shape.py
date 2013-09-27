import numpy as np
import getSummit as gs
import wig as WIG
import argparse
from operator import itemgetter
from sortedcollection import SortedCollection
import bisect
from dpgmm import dpgmm

class Peak:
    '''
    A class to represent peaks with various functions such as dissect peaks etc.
    '''
    def __init__(self, ptuple):
        self.chrom, self.start, self.end, self.name, self.score, self.strand = ptuple
        self.dpgmmModel = None
        self.dpgmmDraw = None
        self.subPeaks = None
        self.totalTags = None
        self.wig = None

    def setWig(self, chromWig):
        '''
        Set the original wig data within the peak region.
        Argument: chromWig, the wig data for the whole chromosome.
        Result: the self.wig will be set with the same columns as the chromWig.
        '''
        ws = bisect.bisect_left(chromWig[:,0], self.start)
        we = bisect.bisect_right(chromWig[:,0], self.end)
        self.wig = chromWig[ws:we,:]
        self.totalTags = self.wig[:,1].sum()

    def __add_data_to_model(self):
        '''
        Helper method to populate the DPGMM model with the wig data.
        '''
        for i in range(self.wig.shape[0]):
            for j in range(int(round(self.wig[i,1]))):
                self.dpgmmModel.add([ self.wig[i,0] ])

    def dissect(self, n_components=5):
        '''
        Apply DPGMM to the peak and find the components.
        Argument: n_components (optional), the maximum number of components. Defaults to 5, which is sufficient in most cases.
        '''
        self.dpgmmModel = dpgmm.DPGMM(1, n_components)
        self.__add_data_to_model()
        self.dpgmmModel.setPrior()
        self.dpgmmModel.solve()
        self.dpgmmDraw = np.array(self.dpgmmModel.intMixture())

    def selectComponents(self, frac=0.3):
        '''
        Returns an array of the bool values of whether a component is selected.
        The threshold is the max_weight * frac.
        Argument: frac, the fraction of the maximum weight to be used as the threshold. Defaults to 0.3.
        '''
        return self.dpgmmDraw[0] > self.dpgmmDraw[0].max() * frac


    def populateSubPeaks(self):
        '''
        Calculate the subpeaks given the mixture model.
        '''
        sbools = self.selectComponents()
        sweights = self.dpgmmDraw[0][sbools]

        scaleFactor = 1 / sweights.sum()

        selected = zip(scaleFactor * self.dpgmmDraw[0][sbools], self.dpgmmDraw[1][sbools])
        selected.sort(key = lambda k: k[1].getLoc() )
        maxs, mins = self._getLocalMaxAndMin(selected)
        lastMin = 0
        for m in mins:
            self.subPeaks.append(Peak((self.chrom, self.start + lastMin, self.start + m, self.name+'_1', self.score, self.strand  )))
            lastMin = m

    def _getProb(self, selected,x):
        '''
        Return the probability density of a point.
        '''
        y = 0
        for s in selected:
            y += s[0] * s[1].prob(x)

        return x

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
        maxs = argrelextrema(y, np.greater)[0]
        mins = argrelextrema(y, np.less)[0]

        return maxs, mins










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
            if strand == '-':
                if chrom not in rpeaks:
                    rpeaks[ chrom ] = []
                rpeaks[chrom].append( [tokens[0], start, end, name, score, strand] )
            else:
                if chrom not in fpeaks:
                    fpeaks[ chrom ] = []
                fpeaks[chrom].append( [tokens[0], start, end, name, score, strand] )
    f.close()

    for chrom in fpeaks:
        fpeaks[chrom].sort(key = lambda k:( k[0], k[1], k[2] ))
        if chrom in rpeaks:
            rpeaks[chrom].sort(key = lambda k:( k[0], k[1], k[2] ))
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



def pair( fpeaks, rpeaks, fwig, rwig, ulimit, dlimit, prefix):
    '''
    Assuming that the peaks on one strand is mutually exclusive.
    They do not overlap with each other.
    In this case, the ordering of the starts of the peaks and the
    ends of the peaks are the same. And that when the starts are
    sorted, the ends are also sorted.
    '''
    #print "fwig: ", fwig
    #print "rwig: ", rwig
    offset = 5
    expandCol = 1
    out1 = open(prefix + "_singletons_shape_c2c.bed",'w')
    out2 = open(prefix + "_pairs_shape_c2c.gff", "w") #c2c means center-to-center
    out3 = open(prefix + "_pairs_shape_true_e2e_detail.gff",'w') #e2e means end-to-end.
    out4 = open(prefix + "_pairs_shape_true_e2e.gff",'w')
    print fpeaks.keys()
    for chrom in fpeaks:
        if chrom not in rpeaks:
            continue
        print chrom
        fp = fpeaks[chrom]
        rp = rpeaks[chrom]
        pairF = []  #Store the pairing information, if unpaired, it will be negative.
        pairR = []
        fw = fwig[chrom]
        expandedFw = WIG.expandWig( fw, offset, expandCol, smooth=False )
        rw = rwig[chrom]
        expandedRw = WIG.expandWig( rw, offset, expandCol, smooth=False )
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
            rstarts.append( r[1] )
            rends.append( r[2] )
        for i in range( len( fp ) ):
            currfp = fp[ i ]
            start = currfp[1]
            end = currfp[2]
            es = start - ulimit
            ee = end + dlimit
            currFw = expandedFw[ max( 0, start - fw[0,0] ) + offset : max( 0, end - fw[0,0] ) + offset + 1]
            si = bisect.bisect_left( rends, es )
            ei = bisect.bisect_right( rstarts, ee )
            ftagCounts,_,_ = gs.getTagCount( fwig, chrom, start, end )

            #print ei - si
            maxScore = 0
            bestDist = 0
            bestIdx = 0
            bestRpos = 0
            bestAve = 0
            bestHeight = 0
            for idx in range( si, ei ):
                currrp = rp[ idx ]
                rstart = currrp[1]
                rend = currrp[2]
                currRw = expandedRw[ max( 0, rstart - rw[0, 0] ) + offset : max( 0, rend - rw[0,0] ) + offset + 1 ]
                rtagCoungs,_,_ = gs.getTagCount( rwig, chrom, currrp[1], currrp[2] )

                tempScore, tempDist, tempRpos, tempAve, tempHeight = getScore( currFw, currRw, start, rstart )

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
            fp[i][1] -= 1
            if f[0] == -1:
                singletons.append( fp[i] )
            else:
                rp[f[0]][1] -= 1
                #pairs.append( fp[i] )
                #pairs.append( rp[f[0]] )
                #pairStart = (2*f[3] - f[2])/2
                #pairEnd = pairStart + 1
                pairStart = f[3] - f[2] + 1
                pairEnd = f[3]
                pairs.append( [fp[i][0],f[5],'.',pairStart, pairEnd,f[4],'.','.','cw_distance='+str(f[2]) ] )
                pairs_extended.append( [fp[i][0],'.','.',fp[i][1]+1, rp[f[0]][2],f[5],'.','.','bestAve=%f,maxScore=%f,bestDist=%d'%(f[4],f[1],f[2]) ] + fp[i] + rp[f[0]] )

        for i,f in enumerate(pairR):
            rp[i][1] -= 1
            if f[0] == -1:
                singletons.append( rp[i] )

        singletons.sort(key=lambda k:( k[0], k[1], k[2] ))
        pairs.sort(key = lambda k:( k[0], k[1], k[2]))
        print "singletons: ", len(singletons)
        print "pairs: ", len(pairs)

        for s in singletons:
            out1.write('\t'.join([str(i) for i in s]))
            out1.write('\n')
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
    out1.close()
    out2.close()
    out3.close()
    out4.close()



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

    pair( fpeaks, rpeaks, fwig, rwig, args.upstream, args.downstream, args.output)

if __name__ == "__main__":
    main()
