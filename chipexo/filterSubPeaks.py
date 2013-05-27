
#Standard libraries
import os,sys
import argparse
import bisect
from multiprocessing import Process, Queue, current_process, freeze_support

#Third party libraries
import numpy as np

import wig
import pair_shape as ps

def getPeaks( filename ):
    peaka, peakb = ps.loadPeaks( filename )
    for c in peakb:
        if c in peaka:
            peaka[c] += peakb[c]
        else:
            peaka[c] = peakb[c]

    for c in peaka:
        peaka[c].sort(key=lambda k:( k[1], k[2] ) )

    return peaka

def filter_a_chrom(chrom, oriPeaks, subPeaks, chromWig, args):
    '''
    Filters a chromsome
    oriPeaks and subPeaks are lists of peaks from the same chromosome.
    Both peak lists are sorted by coordinates on the chromosome.
    '''
    sub_peak_starts = [k[1] for k in subPeaks]
    result = []
    for op in oriPeaks:
        startIdx = bisect.bisect_left( chromWig[:,0], op[1] )
        endIdx = bisect.bisect_right( chromWig[:,0], op[2] )
        subStartIdx = bisect.bisect_left( sub_peak_starts, op[1] )
        subEndIdx = bisect.bisect_right( sub_peak_starts, op[2] )
        if endIdx <= startIdx or subEndIdx <= subStartIdx:
            continue
        currWig = wig.expandWig( chromWig[startIdx:endIdx,:],0, 1, smooth=False )
        startPos = chromWig[ startIdx, 0 ]
        endPos = chromWig[ endIdx - 1, 0 ]
        maxV = np.max( currWig ) #The maximum value
        for subI in range( subStartIdx, subEndIdx ):
            relStart = max( 0, subPeaks[ subI ][1] - startPos )
            relEnd = min( currWig.shape[0] - 1, subPeaks[ subI ][2] - startPos )
            subMaxV = np.max( currWig[ relStart: relEnd ] )
            if subMaxV > args.threshfrac * maxV and subMaxV > args.cutoff:
                result.append( subPeaks[ subI ] )
                result[-1][4] = subMaxV
                result[-1][5] = args.strand
    return result

def filter_a_chrom_process( tasks, args, outQ ):
    '''
    An intermediate method for processing the tasks queue and call filter_a_chrom
    '''
    print 'Filtering subpeaks.'
    for chrom, oriPeaks, subPeaks, chromWig in iter(tasks.get, "STOP"):
        print 'Process ', processID, ' is processing ', chrom
        outQ.put( (chrom, filter_a_chrom( chrom, oriPeaks, subPeaks, chromWig, args ) ) )
        

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Filter splitted peaks reported by PeakSplitter.')
    parser.add_argument('wigfile', help='Wig file for the peaks.')
    parser.add_argument('subpeakfile', help='Peak file in bed format.')
    parser.add_argument('oripeakfile', help='The original bed file contains the peaks splitted.' )
    parser.add_argument('strand', choices=['+', '-', '.'], help='The strand of the peaks we are currently processing. Choices are +, -, . .')
    parser.add_argument( '--threshfrac', type=float, default = 0.5, help='minimum pileup of a subpeak as a fraction of the highest point of the original peak. Default is 0.5.' )
    parser.add_argument('--cutoff', type=float, default = 5, help='minimum pileup of a subpeak regardless the relative fraction. Default is 5.')

    args = parser.parse_args()
    wigData = wig.loadWig( args.wigfile, smooth = False )
    subPeaks = getPeaks( args.subpeakfile )
    oriPeaks = getPeaks( args.oripeakfile )
    chroms = set( subPeaks.keys() ).intersection(set( oriPeaks.keys() ))
    tasks = Queue()
    count = 0
    for c in chroms:
        count += 1
        tasks.put( (c, oriPeaks[c], subPeaks[c], wigData[c]) )

    freeze_support()
    NUM_PROCESSES = 4
    processID = 1
    processes = []
    outQ = Queue()
    for i in range( NUM_PROCESSES ):
        processes.append( Process( target = filter_a_chrom_process, args=( tasks, args, outQ )))
        processes[-1].start()
        processID += 1
    outTokens = args.subpeakfile.split('.')
    outTokens[-1] = 'filtered.' + outTokens[-1]
    out = open('.'.join(outTokens), 'w')
    for i in range( count ):
        (chrom, resultPeaks) = outQ.get()
        print 'writing ', chrom
        for p in resultPeaks:
            out.write('\t'.join([str(k) for k in p]))
            out.write('\n')

    for i in range( NUM_PROCESSES ):
        tasks.put('STOP')
    for i in range( NUM_PROCESSES ):
        processes[i].join()

    out.close()
