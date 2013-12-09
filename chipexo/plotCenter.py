import numpy as np
import pylab as plt
import argparse
import getSummit as gs

def loadPos( filename, chromCol, startCol, endCol, strandCol, offset, coordFormat ):
    f = open( filename )
    pos = {}

    lastChrom = ''
    lastPos = 0
    lastStrand = ''
    for r in f:
        tokens = r.strip().split('\t')
        chrom = tokens[ chromCol ].strip()
        start = int( tokens[ startCol ].strip() )
        end = int( tokens[ endCol ].strip() )
        strand = tokens[ strandCol ].strip()
        currPos = start + offset
        if coordFormat == 0:
            currPos += 1
        if strand == '-':
            currPos = end - offset
        if chrom not in pos:
            pos[ chrom ] = set()
        if lastChrom != chrom or lastPos != currPos or lastStrand != strand:
            pos[ chrom ].add( (currPos,strand) )
            lastChrom = chrom
            lastPos = currPos
            lastStrand = strand

    for chrom in pos:
        pos[chrom] = list(pos[chrom])
        pos[chrom].sort(key=lambda k:(k[0], k[1],))

    f.close()
    return pos

def writeAll( values, outName ):
    f = open( outName, 'w')
    for r in values:
        f.write( '\t'.join([str(i) for i in r ] ))
        f.write('\n')
    f.close()


def plot( fvalues, rvalues, width, out ):
    x = range(-width, width+1)
    for f in fvalues:
        plt.plot(x, f, '--', color='blue')
    for r in rvalues:
        plt.plot(x, r, '--', color='red')
    plt.savefig(out+".png", dpi=600)
    #plt.show()

def getMappedCount(fwig, rwig):
    mappedCount = 0
    for chrom in fwig:
        mappedCount += fwig[chrom][:,1].sum()
    for chrom in rwig:
        mappedCount += rwig[chrom][:,1].sum()
    return mappedCount

def main(args):
    fwig = gs.loadWig( args.forwardWig, smooth=False )
    rwig = gs.loadWig( args.reverseWig, smooth=False )
    mappedCount = getMappedCount(fwig, rwig)
    poses = loadPos( args.positionFile, args.chromCol, args.startCol, args.endCol, args.strandCol, args.offset, args.format )
    values = []
    #print "\n"
    #print fwig.keys()
    #print rwig.keys()
    for chrom in poses:
        #print "\n"
        #print chrom
        if chrom in fwig:
            chromFwig = gs.expandWig( fwig[ chrom ], 0, 1, False )
        if chrom in rwig:
            chromRwig = gs.expandWig( rwig[ chrom ], 0, 1, False)
        chromPos = poses[ chrom ]
        for p,strand in chromPos:
            keep = True
            tempFValues = np.zeros( 2 * args.width + 1 )
            tempRValues = np.zeros( 2 * args.width + 1 )
            if chrom in fwig:
                #print chromFwig.shape
                start = int(p - args.width - fwig[ chrom ][0,0])
                end = int(p + args.width - fwig[ chrom ][0,0])
                #print 'Forward'
                #print start, " " , end
                #print abs(min(0,start)), ' ', min( tempFValues.shape[0], tempFValues.shape[0] + chromFwig.shape[0] - end -1 )
                #print max(0, start) , ' ', min( chromFwig.shape[0], end + 1 )
                if end >= 0 and start < chromFwig.shape[0]:
                    tempFValues[abs(min(0,start)):min( tempFValues.shape[0], tempFValues.shape[0] + chromFwig.shape[0] - end -1) ] = chromFwig[ max(0, start) : min( chromFwig.shape[0], end + 1 ) ]
                    if tempFValues.sum() < args.thresh:
                        keep = False

            if chrom in rwig:
                #print chromRwig.shape
                start = int(p - args.width - rwig[ chrom ][0,0])
                end = int(p + args.width - rwig[ chrom ][0,0])
                #print 'Reverse'
                #print start, " ", end
                #print abs(min(0,start)), ' ', min( tempRValues.shape[0], tempRValues.shape[0] + chromRwig.shape[0] - end -1)
                #print max(0, start) , ' ', min( chromRwig.shape[0], end + 1 )
                if end >= 0 and start < chromRwig.shape[0]:
                    tempRValues[abs(min(0,start)):min( tempRValues.shape[0], tempRValues.shape[0] + chromRwig.shape[0] - end -1 ) ] = chromRwig[ max(0, start) : min( chromRwig.shape[0], end + 1 ) ]
                    if tempRValues.sum() < args.thresh:
                        keep = False
            if keep:
                if strand == '-':
                    temp = tempFValues[::-1]
                    tempFValues = tempRValues[::-1]
                    tempRValues = temp
                values.append(10**6*np.array(np.ma.concatenate([tempFValues, tempRValues])) / mappedCount)


    values.sort(key=lambda k:(sum(k),))

    writeAll( values, args.out + ".txt" )
    values = np.array(values)
    plot( values[:,0:2*args.width+1], values[:,2*args.width+1:], args.width, args.out )



if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Plot tag distribution around the centered positions")
    parser.add_argument("forwardWig", help="The wig file for the forward strand")
    parser.add_argument("reverseWig", help="The wig file for the reverse strand")
    parser.add_argument("positionFile", help="The file that contains the information")
    parser.add_argument("chromCol", type=int, help="The column that contains the chrom information. The column's index start from 0.")
    parser.add_argument("startCol", type=int, help="The column that contains the start positions. The column's index start from 0.")
    parser.add_argument("endCol", type=int, help="The column that contains the end positions. The column's index start from 0.")
    parser.add_argument("strandCol", type=int, help="The column that contains the strand information. The column's index start from 0.")
    #parser.add_argument("mappedCounts", type=int, help="The count of mapped reads for calculation of rpm.")
    parser.add_argument('--out', default="", help="The prefix for the output.")
    parser.add_argument('--width', type=int, default=100, help="How much on each side of the positions to plot. Default: 100.")
    parser.add_argument('--offset', type=int, default=0, help="The offset from the start position. Default: 0.")
    parser.add_argument('--format', type=int, default=0, help="The coordinate format. It can be 0, 1.\n0: 0-based, half-open (end is 1-based).\n1: 1-based, closed interval.\nDefault: 0.")
    parser.add_argument('--thresh', type=float, default=0, help="The minimum number of tags required to write to matrix.")

    args = parser.parse_args()
    main(args)
