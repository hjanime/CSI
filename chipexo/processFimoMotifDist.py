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
    parser.add_argument('-p', '--plotsdir', help="Directory to save the plots.")

    return parser.parse_args()

def getRecords( infile ):
    '''
    Input is a GFF file of fimo output.
    Returns a non-dup list of records, a list of motifs, a dict of motif distribution.
    '''
    f = open( infile )
    records = []
    motifs = set()
    details_poses = {}
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
        tempRecord = [seqname.upper(), attrs, dist, motifName.upper()]
        motifs.add( tempRecord[-1] )
        pos_rec = [ int(tokens[3]), int(tokens[4]), tokens[6] ]

        if tempRecord[-1] not in details_poses:
            #details_poses[ tempRecord[-1] ] = { tempRecord[0]:[ pos_rec,],}
            details_poses[ tempRecord[-1] ] = {}
        if tempRecord[0] not in details_poses[ tempRecord[-1] ]:
            details_poses[ tempRecord[-1] ][ tempRecord[0] ] = []

        details_poses[ tempRecord[-1] ][ tempRecord[0] ].append( pos_rec )
        records.append( tempRecord )
    f.close()
    motifs = list( motifs )
    motifs.sort()
    records.sort(key=lambda k:(k[-1],k[2],k[0]))
    records_nodup = []
    assert len(records) > 0
    records_nodup.append( records[0] )
    for r in records[1:]:
        if records_nodup[-1][0] == r[0] and records_nodup[-1][-1] == r[-1]:
            continue
        else:
            records_nodup.append(r)

    for m in motifs:
        for sn in details_poses[ m ]:
            details_poses[ m ][ sn ].sort(key=lambda k:(k[0],k[1],))
            temp_detail_nodup = [details_poses[m][sn][0],]
            for d in details_poses[ m ][ sn ][1:]:
                if d[0] > temp_detail_nodup[-1][1]:
                    temp_detail_nodup.append( d )
            details_poses[ m ][ sn ] = temp_detail_nodup

    return records_nodup, motifs, details_poses

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

def printDistBox( distBtMotifs, intervals, motifs, translate, plotDir ):
    from plotDistCat import plotBox
    for i,mdists in enumerate(distBtMotifs):
        m = motifs[i]
        if translate != None:
            m = translate[ motifs[i] ]
        plotBox( m.replace('.','_'), mdists, intervals, '.', plotDir)


def getMotifCenterPos( start, end ):
    '''
    The start and end are 1-based, inclusive.
    '''
    return (start + end)/2


def printLocHist(leftPoses, rightPoses,singlePoses, motifs, translate, savefig = None):
    '''
    Generate a plot of the distribution of the left and right motifs within each sequence.
    '''
    import pylab as pl
    import numpy as np
    from plotDistCat import getHistBinCenters, saveOrPrint
    for i,m in enumerate( motifs ):
        temp = [ (a,b) for a,b in zip(leftPoses[i], rightPoses[i])]
        temp.sort(key=lambda k:(k[0],k[1]))
        fig = pl.figure()
        subp = fig.add_subplot(211)
        if translate != None and m in translate:
            m = translate[m]
        subp.set_title(m)
        subp.plot([k[0] for k in temp], [k[1] for k in temp],'.')
        ly, binCenters = getHistBinCenters( leftPoses[i], 40, normed=1)
        lr, binCenters = getHistBinCenters( rightPoses[i], 40, normed=1)
        ls, binCenters = getHistBinCenters( singlePoses[i], 40, normed=1)
        subp.set_xlabel("Left position")
        subp.set_ylabel("Right position")
        subp2 = fig.add_subplot(212)
        subp2.plot( binCenters, ly, '-', label="Left (%d)"%len(leftPoses[i]))
        subp2.plot( binCenters, lr, '-', label="Right (%d)"%len(rightPoses[i]))
        subp2.plot( binCenters, ls, '-', label="Single (%d)"%len(singlePoses[i]))
        subp2.legend()

        subp2.set_xlabel("Position")
        subp2.set_ylabel("Frequency")
        saveOrPrint(fig, "LocHist_"+ m, '', savefig)
        pl.close()

def generateHeatMatrix( motifs, details, records, seq_length, intervals , savefig=None):
    import numpy as np
    import pylab as pl
    from plotDistCat import saveOrPrint, getIntervalIdx
    data = {}
    motif2value = {}
    for i,m in enumerate(motifs):
        print m
        motif2value[ m ] = i

    for r in records:
        if r[-1] not in motif2value:
            continue
        if r[0] not in data:
            data[r[0]] = [r[2],]
            for i in range(len(motifs)):
                data[r[0]].append([])
        for p in details[r[-1]][r[0]]:
            data[r[0]][motif2value[r[-1]]+1].append(( p[0],p[1]))

    matrix_all = []
    matrix_exist = []

    count_diff = 0
    count_all = 0
    count_all_by_dist = []
    count_exist_by_dist = []
    #the following are used to store the indexes of the rows in 
    # the complete set
    matrix_by_dist = []
    matrix_exist_by_dist = []
    for i in range(len(intervals) - 1):
        matrix_by_dist.append([])
        matrix_exist_by_dist.append([])
        count_all_by_dist.append(0)
        count_exist_by_dist.append(0)
    NUMPLOTS = len(intervals)


    for seq in data:
        exist = True
        temp = np.zeros(seq_length) - 1
        for i,p in enumerate(data[seq][1:]):
            for t in p:
                temp[t[0]-1:t[1]] = i
            if len(p) == 0:
                exist = False
        temp = list(temp)
        if sum(temp) == -1*seq_length:
            continue
        interIdx = getIntervalIdx(intervals, data[seq][0])

        matrix_all.append(temp)
        if interIdx != None:
            matrix_by_dist[interIdx].append(len(matrix_all)-1)
            count_all_by_dist[interIdx] += 1

        count_all += 1
        if exist:
            matrix_exist.append(temp)
            if interIdx != None:
                matrix_exist_by_dist[interIdx].append(len(matrix_exist)-1)
                count_exist_by_dist[interIdx] += 1
        else:
            count_diff += 1



    print count_diff, ' ', count_all, ' ', len(matrix_all),' ', len(matrix_exist)
    print count_all_by_dist
    print count_exist_by_dist
    fig_all = pl.figure(1)
    fig_exist = pl.figure(2)
    sub_count = 1
    matrix_all_idx = range(len(matrix_all))
    matrix_exist_idx = range(len(matrix_exist))
    matrix_all_idx.sort(key=lambda k:(np.average(matrix_all[k]),))
    matrix_exist_idx.sort(key=lambda k:(np.average(matrix_exist[k]),))
    for i in range(len(intervals) - 1):
        matrix_by_dist[i].sort(key=lambda k:(np.average(matrix_all[k]),))
        matrix_exist_by_dist[i].sort(key=lambda k:(np.average(matrix_exist[k]),))


    subTitles = []
    for i in range(len(intervals) - 1):
        t = '[%d, %d'%(intervals[i], intervals[i+1])
        if i != len(intervals) - 2:
            t += ')'
        else:
            t += ']'
        subTitles.append(t)

    matrix_all = np.array(matrix_all)
    matrix_exist = np.array(matrix_exist)
    matrix_all[matrix_all<0] = None
    matrix_exist[matrix_exist<0] = None
    for i in range(len(intervals) - 1):
        fas = fig_all.add_subplot((NUMPLOTS+1)/2, 2, sub_count)
        fes = fig_exist.add_subplot((NUMPLOTS+1)/2, 2, sub_count)
        fas.imshow(matrix_all[matrix_by_dist[i]], aspect='auto')
        fas.set_title(subTitles[ sub_count - 1])
        fes.imshow(matrix_exist[matrix_exist_by_dist[i]], aspect='auto')
        fes.set_title(subTitles[ sub_count - 1])
        sub_count += 1

    #matrix_all.sort(key=lambda k:(np.average(k),)) # np.average(k[0:60]), np.average(k[60:140]),np.average(k[140:])))
    #matrix_exist.sort(key=lambda k:(np.average(k),)) # np.average(k[0:60]), np.average(k[60:140]),np.average(k[140:])))

    fas = fig_all.add_subplot((NUMPLOTS+1)/2, 2, sub_count)
    fes = fig_exist.add_subplot((NUMPLOTS+1)/2, 2, sub_count)
    fas.imshow(matrix_all[matrix_all_idx], aspect='auto')
    fas.set_title("All dist")
    fig_all.suptitle("Heatmap all for CTCF (blue) and MAX (red)")
    #pl.colorbar()
    saveOrPrint(fig_all, "Heatmap_all", '', savefig)
    #fig_all.close()
    pl.close()

    fes.imshow(matrix_exist[matrix_exist_idx], aspect='auto')
    fig_exist.suptitle("Heatmap exist for CTCF (blue) and MAX (red)")
    #pl.colorbar()
    fes.set_title("All dist")
    saveOrPrint(fig_exist, "Heatmap_exist", '', savefig)
    #fig_exist.close()
    pl.close()



def main( args ):
    records, motifs, details_poses = getRecords( args.input )
    args.intervals = list(set(args.intervals))
    args.intervals.sort()
    if len( args.intervals ) < 2:
        print "Intervals not defined."
        return
    result = []
    distBtMotifs = []
    leftPoses = []
    rightPoses = []
    singlePoses = []
    for i in range( len(motifs) ):
        temp = []
        distTemp = []
        for j in range( len(args.intervals) - 1 ):
            temp.append({})
            distTemp.append([])
        result.append(temp)
        distBtMotifs.append( distTemp )
        leftPoses.append([])
        rightPoses.append([])
        singlePoses.append([])

    if len(records) < 1 : sys.exit()
    motifIndex = 0
    currMotif = motifs[motifIndex]
    currSeqName = records[0][0]
    currSeqMotifCount = 0
    intervalIndex = getIntervalIdx( args.intervals, records[0][2])
    print args.intervals
    for r in records:
        detail_r = details_poses[r[-1]][r[0]]
        if len(detail_r) > 1:
            if intervalIndex != None:
                tempDist = detail_r[-1][0] - detail_r[0][0]
                distBtMotifs[ motifIndex ] [intervalIndex].append( tempDist )
            leftPoses[ motifIndex ].append( getMotifCenterPos( detail_r[0][0], detail_r[0][1] ) )
            rightPoses[ motifIndex ].append( getMotifCenterPos( detail_r[-1][0], detail_r[-1][1] ) )
        else:
            singlePoses[ motifIndex ].append( getMotifCenterPos( detail_r[0][0], detail_r[0][1] ) )
        if r[-1] != currMotif:
            motifIndex += 1
            assert r[-1] == motifs[motifIndex]
            currMotif = motifs[motifIndex]
        else:
            currSeqMotifCount = len( details_poses[r[-1]][r[0]] )
            #print r[2],' ', intervalIndex
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


    printDistBox( distBtMotifs, args.intervals, motifs, translate, args.plotsdir)
    printLocHist( leftPoses, rightPoses,singlePoses, motifs, translate, args.plotsdir )

    generateHeatMatrix( ['MA0139.1','MA0058.1'], details_poses, records, 201, args.intervals, args.plotsdir)

if __name__=='__main__':
    args = parseArgs()
    main(args)
