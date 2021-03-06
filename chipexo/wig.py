import os,sys
import numpy as np
import scipy.stats as scist
import scipy.misc as scimi
import gc

from multiprocessing import Process, Queue, current_process, freeze_support


def density( expanded, sd, nsd, oneOnly=False ):
    '''
    expanded: the 1-dimensional array with all the values of a chromosome
    sd: standard deviation
    nsd: the number of standard deviations to use as cutoff
    '''
    out = np.zeros_like( expanded, dtype=np.float64 )
    for i in range( expanded.shape[0] ):
        count = expanded[i]
        if oneOnly and count > 1:
            count = 1
        if count > 0:
            window = sd * nsd * count
            start = max( 0, i - window )
            end = min( expanded.shape[0], i + window )
            for j in range(int( start ), int( end ) ):
                beta = ( j - i ) * 1.0 / sd
                out[j] += count * 1.0 * np.exp( -0.5 * beta * beta )
    expanded.resize(100000, refcheck=False)
    expanded.resize(0, refcheck=False)
    return out

def density_nb( expanded, r, mean, strand ):
    '''
    expanded is the list of values.
    r, p follows the description in http://en.wikipedia.org/wiki/Negative_binomial_distribution
    mean = p*r/(1-p)
    mode = floor( p(r-1)/(1-p) )
    '''
    mean = float( mean )
    r = float(r)
    p = mean / (mean + r)
    mode = np.floor( p*(r-1)/(1-p ) )
    p = 1 - p #conform to the scipy definition
    nbinom = scist.nbinom( r, p )
    modep = nbinom.pmf( mode )
    factor = 1/modep
    leftwin = mode
    rightwin = mean * 3
    if strand == '-':
        temp = leftwin
        leftwin = rightwin
        rightwin = temp

    out = np.zeros_like( expanded )
    for i in range( expanded.shape[0] ):
        count = expanded[i]
        if count > 0:
            start = max( 0, i - leftwin )
            end = min( expanded.shape[0], i + rightwin )
            for j in range( int(start), int(end) ):
                k = j - i + mode
                if strand == '-':
                    k = mode - j + i
                out[ j ] += factor*count*nbinom.pmf( k )
    expanded.resize( 100000, refcheck=False)
    expanded.resize( 0, refcheck=False)

    return out




def expandWig( chromWig, offset, expandCol, smooth=True, strand = '+', method='g', kargs=None ):
    '''
    strand is not used in Gaussian smoothing.
    '''
    assert chromWig != None
    assert chromWig.shape[0] > 0 and chromWig.shape[1] > 0 and expandCol < chromWig.shape[1]
    assert offset >= 0
    assert method=='g' or method =='nb' or method=='go'
    assert kargs == None or len( kargs ) > 1
    if not kargs:
        if method == 'g':
            kargs = (3,3) # bandwidth, number of bandwidth
        elif method == 'nb':
            kargs = (2, 10) #r, mean
    startp = chromWig[ 0, 0 ]
    endp = chromWig[ -1, 0 ]
    expanded = np.zeros( endp - startp + 1 + 2*offset )
    for i in range( chromWig.shape[0] ):
        expanded[ chromWig[i, 0] - startp + offset ] = chromWig[ i, expandCol ]

    if smooth:
        if method == 'g':
            expanded = density( expanded, kargs[0], kargs[1])
        elif method == 'go':
            expanded = density( expanded, kargs[0], kargs[1], True)
        else:
            expanded = density_nb( expanded, kargs[0], kargs[1], strand)
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
    lastChrom = None
    lastSpan = -1
    currSpan = 1
    for r in f:
        if r.startswith('variable'):
            lastChrom = currChrom
            lastSpan = currSpan

            tokens = r.strip().split()
            currChrom = tokens[1].split('=')[1]
            currSpan = 1
            if len(tokens) > 2 and tokens[2].startswith("span"):
                currSpan = int(tokens[2].split('=')[1])

            if currChrom == lastChrom and currSpan == lastSpan:
                continue
            if len( temp ) > 0 and lastChrom != None :
                lines.put( (lastChrom, temp, lastSpan) )
                count += 1
            temp = []
        elif r.startswith('track'):
            continue
        else:
            temp.append( r.strip() )
    if len(temp) > 0:
        lines.put( (currChrom, temp, currSpan) )
        count += 1
    f.close()
    return lines, count

def processChrom( taskQ, outQ, processID, offset, smooth = True, strand = '+' ):
    for chrom,chromLines,span in iter(taskQ.get, "STOP"):
        print 'Process ', processID, ' is processing ', chrom, ' span is ', span
        temp = [ ]
        for line in chromLines:
            tokens = line.strip().split()
            #if float(tokens[1]) < 3:
            #    continue
            for i in range(span):
                temp.append([ int(tokens[0]) + i, float( tokens[1] ), 0])
        print 'Process ', processID, ' sorting-------'
        temp.sort( key= lambda k:(k[0]))
        temp = np.array( temp )

        if smooth:
            print "Process ", processID, " is smoothing."
            startp = temp[0, 0]
            forGaussian = expandWig( temp, offset, 1 )
            for i in range( temp.shape[0] ):
                temp[ i, 2 ] = forGaussian[ temp[i,0] - startp + offset ]
            forGaussian.resize(100000, refcheck=False)
            forGaussian.resize(0, refcheck=False)

        outQ.put( ( chrom, temp) )
        #temp = None
        del chromLines
        print 'Process ', processID, ' cleaning'
        gc.collect()


def loadWig(filename, smooth = True, strand = '+'):
    freeze_support()
    lines, count = readWigFile( filename )
    outQ = Queue()

    wig = {}
    print "loading wig --------------"

    NUM_PROCESSES = 3
    if not smooth:
        NUM_PROCESSES = 4
    processID = 1
    processes = []
    offset = 5

    for i in range( NUM_PROCESSES ):
        processes.append(Process( target = processChrom, args = ( lines, outQ, processID, offset, smooth, strand  ) ))
        processes[-1].start()
        processID += 1

    for i in range( count ):
        temp = outQ.get()
        print 'storing ', temp[0], ' ', type(temp[0]), ' ', temp[1].shape
        if temp[0] in wig:
            wig[temp[0]] = np.concatenate((wig[temp[0]], temp[1]), axis=0)
        else:
            wig[ temp[0] ] = temp[1]
    for i in range( NUM_PROCESSES ):
        lines.put( "STOP" )
    for i in range( NUM_PROCESSES ):
        processes[i].join()

    print 'Done loading ', filename
    return wig

def writeAsBedGraph(wig, filename):
    '''
    Write a wig input object as bedGraph format.
    '''
    out = open(filename, 'w')
    chroms = wig.keys()
    chroms.sort()
    for chrom in chroms:
        chromWig = wig[chrom]
        for p in chromWig:
            out.write('%s\t%d\t%d\t%f\n'%(chrom, p[0] - 1 , p[0], p[1]))
    out.close()


if __name__=='__main__':
    wig = loadWig( sys.argv[1], strand = sys.argv[2] )
