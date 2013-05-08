import os,sys
import numpy as np
import gc

from multiprocessing import Process, Queue, current_process, freeze_support


def density( expanded, sd, nsd ):
    '''
    expanded: the 1-dimensional array with all the values of a chromosome
    sd: standard deviation
    nsd: the number of standard deviations to use as cutoff
    '''
    out = np.zeros_like( expanded, dtype=np.float64 )
    for i in range( expanded.shape[0] ):
        count = expanded[i]
        if count > 0:
            window = sd * nsd * count
            start = max( 0, i - window )
            end = min( expanded.shape[0], i + window )
            for j in range(int( start ), int( end ) ):
                beta = ( j - i ) * 1.0 / sd
                out[j] += count * 1.0 * np.exp( -0.5 * beta * beta )
    return out

def expandWig( chromWig, offset, expandCol, smooth=True ):
    assert chromWig != None
    assert chromWig.shape[0] > 0 and chromWig.shape[1] > 0 and expandCol < chromWig.shape[1]
    assert offset >= 0
    startp = chromWig[ 0, 0 ]
    endp = chromWig[ -1, 0 ]
    expanded = np.zeros( endp - startp + 1 + 2*offset )
    for i in range( chromWig.shape[0] ):
        expanded[ chromWig[i, 0] - startp + offset ] = chromWig[ i, expandCol ]

    if smooth:
        expanded = density( expanded, 3, 3)
        #expanded = scf.gaussian_filter1d( expanded,10 )
    gc.collect()
    return expanded


def readWigFile( filename ):
    '''
    Reads a wig file and store the lines in a list by chromosome.
    Each chromosome is a sublist, and the first element is the chromsome name.
    '''
    f = open( filename )
    lines = Queue() 
    count = 0
    temp = []
    currChrom = None
    for r in f:
        if r.startswith('variable'):
            if len( temp ) > 0 and currChrom != None:
                lines.put( (currChrom, temp) )
                count += 1
            temp = []
            currChrom = r.strip().split('=')[1]
        else:
            temp.append( r.strip() )
    if len(temp) > 0:
        lines.put( (currChrom, temp) )
        count += 1
    f.close()
    return lines, count

def processChrom( taskQ, outQ, processID, offset, smooth = True ):
    for chrom,chromLines in iter(taskQ.get, "STOP"):
        print 'Process ', processID, ' is processing ', chrom
        temp = [ ]
        for line in chromLines:
            tokens = line.strip().split()
            temp.append([ int(tokens[0]), float( tokens[1] ), 0])
        print 'Process ', processID, ' sorting-------'
        temp.sort( key= lambda k:(k[0]))
        temp = np.array( temp )
        
        if smooth:
            startp = temp[0, 0]
            forGaussian = expandWig( temp, offset, 1 )
            for i in range( temp.shape[0] ):
                temp[ i, 2 ] = forGaussian[ temp[i,0] - startp + offset ]
            del forGaussian

        outQ.put( ( chrom, temp) )
        #temp = None
        del chromLines
        print 'Process ', processID, ' cleaning'
        gc.collect()


def loadWig(filename, smooth = True):
    freeze_support()
    lines, count = readWigFile( filename )
    outQ = Queue()

    wig = {}
    print "loading wig --------------"

    NUM_PROCESSES = 3
    processID = 1
    processes = []
    offset = 5

    for i in range( NUM_PROCESSES ):
        processes.append(Process( target = processChrom, args = ( lines, outQ, processID, offset, smooth  ) ))
        processes[-1].start()
        processID += 1

    for i in range( count ):
        temp = outQ.get()
        print 'storing ', temp[0]
        wig[ temp[0] ] = temp[1]
    for i in range( NUM_PROCESSES ):
        lines.put( "STOP" )
    for i in range( NUM_PROCESSES ):
        processes[i].join()

    print 'Done loading'
    return wig


if __name__=='__main__':
    wig = loadWig( sys.argv[1] )
