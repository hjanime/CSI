import os,sys
from plotDistCat import getIntervalIdx

def parseArgs():
    import argparse
    parser = argparse.ArgumentParser(description="This script takes a fimo output with the motifs in the last column and distances in the first column to generate a table of number of motifs in different distance categories.")
    parser.add_argument('input', help="The fimo output file")
    parser.add_argument('-I', '--intervals', type=int, nargs='+', help="The edges of distance intervals to split the data into. [) style of edges.")
    parser.add_argument('-o', '--outPrefix', help="The output prefix.")
    parser.add_argument('-f', '--fasta', help="The fasta file that the fimo output is generated from.")
    parser.add_argument('-t', '--translate', help="Translate meme identifier to factor name")

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

def loadTranslate( filename ):
    if not filename:
        return None
    f = open( filename )
    translate = {}
    for r in f:
        if r.startswith('#'):
            continue
        tokens = r.strip().split()
        translate[ tokens[0].upper() ] = tokens[1].upper()

    f.close()
    return translate

def printResult( result, motifs, intervals, outPrefix = None, fastaCount = None, translate=None):
    out = sys.stdout
    outS = None
    if outPrefix:
        out = open( outPrefix + '.detail.tsv','w')
        outS = open( outPrefix + '.tsv','w')
    
    def writeToFile( out, useSum = False ):
        out.write('motifs')
        if translate:
            out.write('\tFactor')
        for i in range(len(intervals) - 1 ):
            if i == len(intervals) - 2:
                out.write('\t[%d,%d]'%(intervals[i], intervals[i+1],))
            else:
                out.write('\t[%d,%d)'%(intervals[i], intervals[i+1],))
        out.write('\n')

        for i in range( len( motifs ) ):
            out.write(motifs[i])
            if translate and motifs[i] in translate:
                out.write('\t' + translate[motifs[i]])
            for j in range( len( intervals )-1 ):
                temp = []
                tempTotal = 0
                for k in result[i][j]:
                    temp.append('%d:%d'%(k,result[i][j][k],))
                    tempTotal += result[i][j][k]
                temp.sort(key=lambda k:(int(k.split(':')[0])))
                if len( temp ) == 0:
                    temp = ['0',]
                if not useSum:
                    out.write('\t%s'%', '.join(temp))
                else:
                    out.write('\t%d'%tempTotal)
            out.write('\n')
        if fastaCount:
            out.write('Number_of_seqs\t')
            for v in fastaCount:
                out.write('\t%d'%v)
            out.write('\n')
    
    writeToFile( out )
    if outPrefix:
        writeToFile( outS, True)
        out.close()
        outS.close()

def getFastaCount( fasta, intervals ):
    f = open( fasta )
    counts = []
    for i in range( len(intervals) -1 ):
        counts.append(0)

    for r in f:
        if r.strip()[0] != '>':
            continue
        v = int(r.split('=')[-1])
        idx = getIntervalIdx( intervals, v)
        if idx != None:
            counts[ idx ] += 1
    return counts

def main( args ):
    records, motifs = getRecords( args.input )
    args.intervals = list(set(args.intervals))
    args.intervals.sort()
    if len( args.intervals ) < 2:
        print "Intervals not defined."
        return
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
                if intervalIndex != None  and currSeqMotifCount in result[motifIndex][intervalIndex]:
                    result[motifIndex][intervalIndex][ currSeqMotifCount ] += 1
                elif intervalIndex != None:
                    result[motifIndex][intervalIndex][ currSeqMotifCount ] = 1
                currSeqMotifCount = 1
                intervalIndex = getIntervalIdx( args.intervals, r[2] )
                currSeqName = r[0]

    fastaCount = getFastaCount( args.fasta, args.intervals )
    translate = loadTranslate( args.translate )
    printResult( result, motifs, args.intervals, args.outPrefix, fastaCount, translate)


if __name__=='__main__':
    args = parseArgs()
    main(args)
