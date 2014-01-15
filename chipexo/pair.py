import os,sys
import numpy as np
import getSummit as gs
import argparse
from operator import itemgetter
from sortedcollection import SortedCollection
import bisect
import scipy.stats as scs


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
            pvalue = -1
            qvalue = -1
            signal = -1
            peak = -1
            if filename.lower().endswith("narrowpeak"):
                signal = float(tokens[6])
                pvalue = float(tokens[7])
                qvalue = float(tokens[8])
                peak = int(tokens[9])

            if strand == '-':
                if chrom not in rpeaks:
                    rpeaks[ chrom ] = []
                rpeaks[chrom].append( [tokens[0], start, end, name, score, strand, signal, pvalue, qvalue, peak] )
            else:
                if chrom not in fpeaks:
                    fpeaks[ chrom ] = []
                fpeaks[chrom].append( [tokens[0], start, end, name, score, strand, signal, pvalue, qvalue, peak] )
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
    fend = fstart + len( fp )
    rend = rstart + len( rp )
    fmean = sum( fp ) / len( fp )
    rmean = sum( rp ) / len( rp )
    newfp = np.array( [ ( t / fmean) ** 2 for t in fp ], dtype="float32")
    newrp = np.array( [ ( t / rmean) ** 2 for t in rp ], dtype="float32")
    maxScore = 0
    shift = 0
    i = 0
    maxRpos = 0
    maxheight = 0
    #Gradually move forward towards the 3' end
    while fstart + i < rend:
        overlap = 0
        tempStart = max( fstart + i, rstart )
        tempEnd = min( fend + i, rend )
        tempRpos = 0
        tempHeight = 0
        for j in range( tempStart, tempEnd ):
            temp =  newfp[ j - ( fstart + i ) ] * newrp[ j - rstart ]
            if temp >  tempHeight:
                tempHeight = temp
                tempRpos = j
            overlap += temp
        #coeff = scs.pearsonr( fp, rp[i:i+len(fp)] )

        tempScore = overlap #coeff #overlap #/ ( abs (sum( rp[ tempStart - rstart : tempEnd - rstart + 1] ) - sum( fp[ tempStart - (fstart+i) : tempEnd - (fstart+i) +1] )) + 10 )
        if tempScore > maxScore:
            maxScore = tempScore
            shift = i
            maxheight = max( fp )
            maxRpos = tempRpos
            #maxRpos = fp.argmax() + i + fstart
        i += 1

    return maxScore, shift, maxRpos



def pair( fpeaks, rpeaks, fwig, rwig, ulimit, dlimit, prefix):
    '''
    Assuming that the peaks on one strand is mutually exclusive.
    They do not overlap with each other.
    In this case, the ordering of the starts of the peaks and the
    ends of the peaks are the same. And that when the starts are
    sorted, the ends are also sorted.
    '''
    offset = 5
    expandCol = 1
    out1 = open(prefix + "_singletons.bed",'w')
    out2 = open(prefix + "_pairs.gff", "w")
    out3 = open(prefix + "_pairs.narrowPeak", "w")
    for chrom in fpeaks:
        if chrom not in rpeaks:
            continue
        print chrom
        fp = fpeaks[chrom]
        rp = rpeaks[chrom]
        pairF = []  #Store the pairing information, if unpaired, it will be negative.
        pairR = []
        fw = fwig[chrom]
        expandedFw = gs.expandWig( fw, offset, expandCol, False )
        rw = rwig[chrom]
        expandedRw = gs.expandWig( rw, offset, expandCol, False )
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
            start = max( start, fw[0,0] )
            end = max( end, fw[0,0] )
            flength = end - start
            si = bisect.bisect_left( rends, es )
            ei = bisect.bisect_right( rstarts, ee )
            ftagCounts,_,_ = gs.getTagCount( fwig, chrom, start, end )
            #print ei - si
            maxScore = 0
            bestDist = 0
            bestIdx = 0
            bestRpos = 0
            for idx in range( si, ei ):
                currrp = rp[ idx ]
                rstart = max(rw[0,0], currrp[1])
                rend = min(rw[-1,0], currrp[2]+flength)
                currRw = expandedRw[ max( 0, start - rw[0, 0] ) + offset : max( 0, rend - rw[0,0] ) + offset + 1 ]
                rtagCoungs,_,_ = gs.getTagCount( rwig, chrom, currrp[1], currrp[2] )

                tempScore, tempDist, tempRpos = getScore( currFw, currRw, start, rstart )

                fprefer[ i ].insert( (idx, tempScore, tempDist, tempRpos) )
                rprefer[ idx ].insert( (i, tempScore, tempDist, tempRpos) )
                if tempScore > maxScore:
                    maxScore = tempScore
                    bestDist = tempDist
                    bestIdx = idx
                    bestRpos = tempRpos
            if maxScore > pairR[ bestIdx ][1]:
                pairF[i] = ( bestIdx, maxScore, bestDist , bestRpos)
                if pairR[ bestIdx ][0] > 0:
                    unpairedF.append( pairR[ bestIdx ][0] )
                pairR[bestIdx] = ( i, maxScore, bestDist, bestRpos )
            else:
                unpairedF.append( i )
            try:
                fprefer[ i ].remove( ( bestIdx, maxScore, bestDist, bestRpos ) )
            except ValueError:
                #print "Value error: ", bestIdx, ' ',maxScore, ' ',bestDist,' ', si,' ', ei
                pass
        singletons = []
        pairs = []
        pairs_narrow = []
        while len(unpairedF) > 0:
            for u in unpairedF:
                if len( fprefer[u] ) > 0:
                    ridx = fprefer[u][-1][0]
                    if pairR[ ridx ][1] < fprefer[u][-1][1]:
                        if pairR[ idx ][0] > 0:
                            pairF[ pairR[ idx ][0] ] = (-1, 0, 0)
                            unpairedF.append( pairR[ idx ][0] )
                        pairR[ ridx ] = ( u, fprefer[u][-1][1], fprefer[u][-1][2], fprefer[u][-1][3] )
                        pairF[ u ] = ( ridx, fprefer[u][-1][1], fprefer[u][-1][2], fprefer[u][-1][3] )
                    fprefer[u].remove( (ridx, fprefer[u][-1][1], fprefer[u][-1][2], fprefer[u][-1][3] ) )
                else:
                    unpairedF.remove( u )

        for i,f in enumerate(pairF):
            fp[i][1] -= 1
            currFp = fp[i]
            if f[0] == -1:
                singletons.append( fp[i] )
            else:
                rp[f[0]][1] -= 1
                currRp = rp[f[0]]
                #pairs.append( fp[i] )
                #pairs.append( rp[f[0]] )
                pairStart = (2*f[3] - f[2])/2
                pairEnd = pairStart + 1
                pairs.append( [fp[i][0],'.','.',pairStart-10, pairEnd+10,f[1],'.','.','cw_distance='+str(f[2]) ] )
                half_len = int(f[2]/2)
                narrowStart=  max(currFp[1] + half_len, currRp[1] - half_len) - half_len
                narrowEnd = min(currFp[2] + half_len, currRp[2] - half_len) + half_len
                pairs_narrow.append([currFp[0], narrowStart, narrowEnd, currFp[3]+'_'+currRp[3], (currFp[4]+currRp[4])/2, '.', (currFp[6]+currRp[6])/2, (currFp[7]+currRp[7])/2, (currFp[8]+currRp[8])/2, f[3] - half_len - narrowStart, f[2]])


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
        for p in pairs_narrow:
            out3.write('\t'.join([str(i) for i in p]))
            out3.write('\n')
    out1.close()
    out2.close()
    out3.close()



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
    fwig = gs.loadWig( args.forwardWig, smooth=False )
    rwig = gs.loadWig( args.reverseWig, smooth=False )

    pair( fpeaks, rpeaks, fwig, rwig, args.upstream, args.downstream, args.output)

if __name__ == "__main__":
    main()
