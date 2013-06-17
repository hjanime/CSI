import os,sys
import pylab as plt
import argparse
import numpy as np
from scipy import stats
from bisect import *
import scipy.ndimage.filters as scf
from wig import *

SMAP = {'+': 'Forward', '-': 'Reverse'}

def loadFragSize(filename):
    f = open(filename)
    fragsizes = {}
    for r in f:
        sample, frag = r.strip().split()
        fragsizes[sample.upper()] = int(frag)
    f.close()
    return fragsizes

def loadReadCount( filename ):
    f = open( filename )
    readCounts = {}
    for r in f:
        sample, count = r.strip().split()
        readCounts[ sample.upper() ] = int( count )
    f.close()
    return readCounts





def getTagCount( wig, chrom, start, end ):

    wigs = bisect_left( wig[chrom][ :, 0 ], start ) # The starting location in the wig
    wige = bisect_right( wig[chrom][:, 0 ], end )
    wigOrig = wig[chrom][ :, 1 ]
    wigv = wig[ chrom ][ : ,2 ]
    tagCount = 0
    numHigh = 0
    summitv = 0
    summitidx = [ wigs ]
    for wi in range( wigs , wige ):
        tagCount += wigOrig[ wi ] #add the number of tags together.
        if wigv[wi] > summitv:
            numHigh = 1
            summitidx[ 0 ] = wi
            summitv = wigv[ wi ]
        elif wigv[ wi ] == summitv:
            if len( summitidx ) <= numHigh:
                summitidx.append( wi )
            else:
                summitidx[ numHigh ] = wi
            numHigh += 1

    return tagCount, numHigh, summitidx, summitv

def loadRejectRegion( filename ):
    f = open( filename )
    reject = {}
    for r in f:
        if not r.startswith('#'):
            tokens = r.strip().split('\t')
            chrom, start, end, name, score, strand = tokens
            start = int(start) + 1 #BED files adjust to 1-based
            end = int(end)
            score = float(score)
            if chrom not in reject:
                reject[ chrom ] = []
            reject[ chrom ].append((chrom, start, end, name, score, strand))
    f.close()
    for key in reject:
        reject[ key ].sort(key = lambda k: (k[1], k[2],))

    return reject

def isReject( reject, chrom, start, end):
    if reject == None:
        return False
    if chrom not in reject:
        return False
    for obj in reject[ chrom ]:
        if min( end, obj[2] ) - max( start, obj[1] ) > 0:
            return True

def writeAll( outxls, peaks, summits, records, thres, filidx, filtered ):
    rejected = 0
    for r in records:
        if float( r[filidx] ) > thres:
            outxls.write( '\t'.join( r ) )
            outxls.write( '\n' )

            peaks.write( '\t'.join( [r[0], str(int(r[1]) - 1), r[2], r[11], r[-2], r[4]] ))
            peaks.write( '\n' )

            for p in r[5].split(','):
                summits.write( '\t'.join( [ r[0], str(int(p)-1), p, r[11], r[-2], r[4] ] ))
                summits.write( '\n' )
        else:
            filtered.write( '\t'.join( [r[0], str(int(r[1]) - 1), r[2], r[11], r[-2], r[4]] ))
            filtered.write( '\n' )

            rejected += 1
    print rejected, " peaks have lower than cutoff rpkm values."


def findTruePeak( args, sample, shift, tsize, strand, outdir, mappedCount, reject=None ):
    commonName = '%s_%s_sh%d_t%d' % ( sample, strand, shift, tsize )
    print commonName
    commonName = os.path.join( outdir, commonName ) # include the outdir in path
    xls = open("%s_peaks.xls"%( commonName, ) )
    
    rejectCount = 0
    count = 0

    outxls = open("%s_peaks.out.xls"%( commonName, ),'w')
    peaks = open("%s_peaks.out.bed"% ( commonName, ), 'w')
    bed = open("%s_summits.out.bed"%( commonName ), 'w')
    filtered = open("%s_peaks.filtered.out.bed"%( commonName, ), 'w')
    readHeader = False
    wig = loadWig(os.path.join(args.wigdir, "%s_%s.wig" % ( sample, SMAP[strand], )), smooth=False)
    records = []
    coverages = []
    lengths = []
    tagCounts = []
    pileup = []
    pileup5 = []
    rpms = []
    for r in xls:
        if r.startswith('#'):
            outxls.write(r)
        elif r.strip() == '':
            outxls.write('#')
            outxls.write(r)
        elif not readHeader:
            outxls.write('#')
            tokens = r.strip().split('\t')
            tokens.insert(4,"strand")
            tokens.insert(7,"5-pileup")
            tokens.append( "tagCount" )
            tokens.append( "RPKM" )
            outxls.write('\t'.join(tokens))
            outxls.write('\n')
            readHeader = True
        else:
            tokens = r.strip().split('\t')
            chrom,start,end,length,summitp,summitv = tokens[0:6]
            start = int(start)
            end = int(end)
            if reject == None or not isReject( reject, chrom, start, end ):
                length = int(length)
                wigs = bisect_left( wig[chrom][ :, 0 ], start ) # The starting location in the wig
                #print 'chromSize ', chrom, ' ', len(wig[chrom]), ' ', wig[chrom][0][-1]
                #print chrom, ' ', wigs, ' ', start
                wige = bisect_right( wig[ chrom ][ :, 0 ], end ) # The end position in the wig, exclusive
                lengths.append( length )
                summitv = float(summitv)
                wigvOrig = wig[ chrom ][ :, 1 ] # the original tag counts.
                wigv = wig[ chrom ][ :, 2 ]
                
                #maxOrig = np.max( wig[ chrom ][ wigs:wige, 1] )
                tagCount, numHigh, summitidx, summitv = getTagCount( wig, chrom, start, end )
                summitps = wig[ chrom ][ :, 0 ][ summitidx[ 0 : numHigh ] ]
                summitp = ','.join( [ str(int(i)) for i in summitps ] )  #get all positions
                tokens[4] = summitp
                pileup.append( float(tokens[5]) )
                pileup5.append( summitv )
                rpkm = tagCount * 1.0 * 1000 * 1000000 / ( length * mappedCount )
                rpm = tagCount * 1.0 * 1000000 / mappedCount
                '''
                for i in range( numHigh ):
                    bed.write( '\t'.join( [ chrom, str(int(summitps[ i ]) - 1), 
                                    str(int(summitps[ i ])), tokens[9], str( summitv ), strand ] ) )
                    bed.write('\n')
                peaks.write( '\t'.join( [ chrom, str(start - 1), 
                                str(end), tokens[9], tokens[8], strand ] ) )
                peaks.write('\n')
                '''
                tokens.insert(4, strand) #add the strand information
                tokens.insert(7, str( summitv )) # The new p
                tokens.append( str( tagCount ) )
                tokens.append( str( rpkm ) )
                #tokens.append( str( rpm ) )
                coverages.append( rpkm )
                rpms.append( rpm )
                tagCounts.append( tagCount )
                if length > 2*shift + 2 and wige - wigs >= 3 and float( tokens[6] ) > min( shift , ( args.pileup * mappedCount * shift ) / 10**7 ): # 2*shift is the pileup frag length, with only one location, it is highly likely to be false
                    records.append(tokens)
                #outxls.write( '\t'.join( tokens ) )
                #outxls.write('\n')
                count += 1
            else:
                rejectCount += 1
    writeAll( outxls, peaks, bed, records, args.rpkm, -1, filtered )
        
    peaks.close()
    bed.close()
    xls.close()
    outxls.close()
    filtered.close()

    print "Peaks changed: ", count
    print "Peaks rejected: ", rejectCount
    return coverages, lengths, tagCounts, pileup, pileup5, rpms


def getShiftAndTsize( fragsize ):
    tsize = 50
    #shift = fragsize / 4 # MACS2 will extend 2*shift during pileup.
    shift = 30
    return tsize, shift

def writeDataHist( data, filename ):
    out = open( filename, 'w' )
    for row in data:
        for x in row:
            out.write(str(x))
            out.write(' ')
        out.write('\n')
    out.close()

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-f', "--fragSize", required=True, help="The file with sample name and fragment sizes.")
    parser.add_argument('-i', "--indir", required=True, help="The directory that has all the sample tags alignments.")
    parser.add_argument("samples", nargs='+', help="Sample names.")
    parser.add_argument('-o', '--outdir', default='.', help="The output directory.")
    parser.add_argument('-w', '--wigdir', required=True, help="The directory for the wig files")
    parser.add_argument('-c', '--callpeak', action='store_true', default=False, help="Whether to call peaks using macs or not")
    parser.add_argument('-b', '--rejectRegion', help="Peaks fall into this region will be discarded. The file must be in BED format with at least 6 fields.")
    parser.add_argument('--rpkm', type=float, default = 0, help="The rpkm threshold used for cutoff of the RPKM values.")
    parser.add_argument('--pileup', type=float, default = 0, help="The pileup threshold to filter the peaks.")
    parser.add_argument('-R', '--readCounts', required=True, help="The file with sample name and read counts that are used for RPKM calculation.")
    parser.add_argument('--shift', type=int, help="The shift size. This will override the values generated by the program." )
    parser.add_argument('--tsize', type=int, help="The tsize size. This will override the values generated by the program." )
    parser.add_argument('--mode', default="auto", help="The keep-dup mode, default is auto. Can be '1', 'all', 'auto'." )
    parser.add_argument('-C','--callOnly',  action="store_true", default=False, help="Call peaks only.")
    args = parser.parse_args()

    fragsizes = loadFragSize(args.fragSize)
    readCounts = loadReadCount( args.readCounts )
    reject = None
    coverages = []
    lengths = []
    tagCounts = []
    labels = []
    pileups = []
    pileup5s = []
    rpms = []
    if args.rejectRegion:
        reject = loadRejectRegion( args.rejectRegion )
    for s in args.samples:
        print "\n\nprocessing ", s
        tsize, shift = getShiftAndTsize( fragsizes[ s.upper() ] )
        if args.shift:
            shift = args.shift
        if args.tsize:
            tsize = args.tsize
        mappedCount = readCounts[ s.upper() ]
        path = os.path.join(args.indir, s)
        if args.callpeak or args.callOnly:
            os.system("macs2 callpeak -B --llocal 20000 --broad --broad-cutoff 0.01 --half-ext -t %s.unique.+.bam -f BAM -g hs -n %s_+_sh%d_t%d --bw 60 --verbose 2 -q 0.001 --nomodel --shiftsize %d --tsize %d --keep-dup %s"%(path, os.path.join( args.outdir, s ), shift, tsize, shift, tsize, args.mode))
            os.system("macs2 callpeak -B --llocal 20000 --broad --broad-cutoff 0.01 --half-ext -t %s.unique.-.bam -f BAM -g hs -n %s_-_sh%d_t%d --bw 60 --verbose 2 -q 0.001 --nomodel --shiftsize %d --tsize %d --keep-dup %s"%(path, os.path.join( args.outdir, s ), shift, tsize, shift, tsize, args.mode))
        if not args.callOnly:
            tempCovs, tempLengths, tempTagCounts, tempPileUp, tempPileUp5, tempRPM = findTruePeak( args, s, shift, tsize, '+', args.outdir, mappedCount, reject ) 
            labels.append( s + "_+" )
            coverages.append( tempCovs )
            lengths.append( tempLengths )
            tagCounts.append( tempTagCounts )
            pileups.append( tempPileUp )
            pileup5s.append( tempPileUp5 )
            rpms.append( tempRPM )
            tempCovs, tempLengths, tempTagCounts, tempPileUp, tempPileUp5, tempRPM = findTruePeak( args, s, shift, tsize, '-', args.outdir, mappedCount, reject )
            labels.append( s + "_-" )
            coverages.append( tempCovs )
            lengths.append( tempLengths )
            tagCounts.append( tempTagCounts )
            pileups.append( tempPileUp )
            pileup5s.append( tempPileUp5 )
            rpms.append( tempRPM )

    if not args.callOnly:

        writeDataHist( coverages, os.path.join( args.outdir, 'coverages_rpkm.txt' ) )

        writeDataHist( lengths, os.path.join( args.outdir, 'lengths.txt' ) )

        writeDataHist( tagCounts, os.path.join( args.outdir, 'tagCounts.txt' ) )

        writeDataHist( pileups, os.path.join( args.outdir, 'pileup.txt' ) )

        writeDataHist( pileup5s, os.path.join( args.outdir, 'pileup5.txt' ) )

        writeDataHist( rpms, os.path.join( args.outdir, "rpm.txt" ) )

        out = open(os.path.join( args.outdir, "labels.txt" ), "w") 
        for l in labels:
            out.write( l )
            out.write( '\n' )
        out.close()

        

        #plt.figure()
        #plt.hist( coverages )
        #plt.show()

if __name__=="__main__":
    main()

