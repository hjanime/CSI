import os,sys

def parseArgs():
    import argparse
    parser = argparse.ArgumentParser(description="This script takes a fimo output with the motifs in the last column and distances in the first column to generate a table of number of motifs in different distance categories.")
    parser.add_argument('input', help="The fimo output file")
    parser.add_argument('-I', '--intervals', type=int, nargs='+', help="The edges of distance intervals to split the data into. [) style of edges.")

    return parser.parse_args()

def getRecords( infile ):
    f = open( infile )
    records = []
    motifs = set()
    for r in f:
        r = r.strip()
        if r[0] == '#':
            continue
        tokens = r.split('\t')
        seqname = tokens[0]
        attrs = tokens[-1]
        dist = int( seqname.split('=')[-1] )
        motifName = attrs.split(';')[0].split('=')[-1]
        motifName = motifName.replace('+','')
        motifName = motifName.replace('-','')
        records.append([seqname.upper(), attrs, dist, motifName.upper()])
        motifs.add( records[-1][-1] )
    f.close()
    motifs = list( motifs )
    motifs.sort()
    records.sort(key=lambda k:(k[-1],k[2],k[0]))
    return records, motifs

def printResult( result, motifs, intervals ):
    sys.stdout.write('motifs')
    for i in range(len(intervals) - 1 ):
        if i == len(intervals) - 2:
            sys.stdout.write('\t[%d,%d]'%(intervals[i], intervals[i+1],))
        else:
            sys.stdout.write('\t[%d,%d)'%(intervals[i], intervals[i+1],))
    sys.stdout.write('\n')

    for i in range( len( motifs ) ):
        sys.stdout.write(motifs[i])
        for j in range( len( intervals )-1 ):
            temp = []
            for k in result[i][j]:
                temp.append('%d:%d'%(k,result[i][j][k],))
            temp.sort(key=lambda k:(int(k.split(':')[0])))
            if len( temp ) == 0:
                temp = ['0',]
            sys.stdout.write('\t%s'%', '.join(temp))
        sys.stdout.write('\n')


def main( args ):
    records, motifs = getRecords( args.input )
    args.intervals = list(set(args.intervals))
    args.intervals.sort()
    if len( args.intervals ) < 2:
        print "Intervals not defined."
        return
    import bisect
    from plotDistCat import getIntervalIdx
    result = []
    for i in range( len(motifs) ):
        temp = []
        for j in range( len(args.intervals) - 1 ):
            temp.append({})
        result.append(temp)

    if len(records) < 1 : sys.exit()
    motifIndex = 0
    currMotif = motifs[motifIndex]
    currSeqName = records[0][0]
    currSeqMotifCount = 0
    intervalIndex = getIntervalIdx( args.intervals, records[0][2])
    for r in records:
        if r[-1] != currMotif:
            motifIndex += 1
            assert r[-1] == motifs[motifIndex]
            currMotif = motifs[motifIndex]
        else:
            if r[0] == currSeqName:
                currSeqMotifCount += 1
            else:
                if intervalIndex  and currSeqMotifCount in result[motifIndex][intervalIndex]:
                    result[motifIndex][intervalIndex][ currSeqMotifCount ] += 1
                elif intervalIndex:
                    result[motifIndex][intervalIndex][ currSeqMotifCount ] = 1
                currSeqMotifCount = 1
                intervalIndex = getIntervalIdx( args.intervals, r[2] )
                currSeqName = r[0]

    printResult( result, motifs, args.intervals )


if __name__=='__main__':
    args = parseArgs()
    main(args)
