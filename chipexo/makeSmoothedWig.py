import os,sys
import wig
import gc

from multiprocessing import Process, Queue, current_process, freeze_support


offset = 0

def add_chrom_data(taskQ, outQ, processID, args, strand  ):
    if args.method == 'g':
        kargs = (args.bw, args.nbw)
    else:
        kargs = (args.r, args.mean)
    for chrom, chromWig in iter( taskQ.get, 'STOP' ):
        print "add chrom data Process ", processID," is processing ", chrom
        lines = []
        lines.append("variableStep chrom=%s\n"%(chrom,))
        startp = chromWig[0,0]
        expanded = wig.expandWig( chromWig, offset, 1, strand = strand, method=args.method, kargs=kargs)
        for i in range( expanded.shape[0] ):
            if expanded[i] > 0.8:
                lines.append( "%d\t%f\n"%(int(i + startp - offset), expanded[i], ) )
        #for i in range( chromWig.shape[0] ):
        #    lines.append('%d\t%f\n'%( chromWig[ i, 0], chromWig[i,1]))
        expanded.resize(100000, refcheck=False)
        expanded.resize(0, refcheck=False)
        chromWig.resize(100000, refcheck=False)
        chromWig.resize(0, refcheck=False)

        outQ.put(lines)
        gc.collect()

def main(args):
    for filename in args.inputs:
        strand = '+'
        if 'Reverse' in filename:
            strand = '-'
        lines = []
        wigdata = wig.loadWig( filename, False, strand )
        print wigdata

        chroms = Queue()
        count = 0
        outQ = Queue()
        for chrom in wigdata:
            print "add ", chrom
            count += 1
            chroms.put( (chrom, wigdata[chrom][:,0:2]) )

        NUM_PROCESSES = args.process

        
        processID = 1
        for i in range( NUM_PROCESSES):
            Process( target=add_chrom_data, args=( chroms, outQ, processID,  args, strand) ).start()
            processID += 1

        if args.method != 'nb':
            kargs = (args.bw, args.nbw)
        else:
            kargs = (args.r, args.mean)
        tokens = filename.split('.')
        tokens[-1] = args.method +'_'+str(kargs[0]) + '_' + str(kargs[1]) +  "_smoothed.wig"
        out = open( '.'.join(tokens), 'w' )
        out.write('track type=wiggle_0 name=%s_%d_%f\n'%( args.method, kargs[0], kargs[1], ))
        for i in range( count ):
            out.write( ''.join(outQ.get()) )


        out.close()
        for i in range( NUM_PROCESSES ):
            chroms.put('STOP')





if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser("Create smoothed wig files")
    parser.add_argument("inputs", metavar="I", nargs='+', help='input wig files')
    parser.add_argument('-m','--method', default='g', choices=['nb', 'g', 'go'], help='kernel for smoothing. nb: negative binomial, g: gaussian, go: gaussian, keep 1 read per position. Default is g')
    parser.add_argument('--bw', type=int, default=3, help='the bandwidth for smoothing, only for Gaussian smoothing')
    parser.add_argument('--nbw', type=int, default=3, help='the number of bandwidths, only for Gaussian smoothing')
    parser.add_argument('--r', type=int, default=2, help='the r in http://en.wikipedia.org/wiki/Negative_binomial_distribution')
    parser.add_argument('--mean', type=float, default=10, help='the mean in http://en.wikipedia.org/wiki/Negative_binomial_distribution')
    parser.add_argument('-p', '--process', type=int, default=3,help='The number of processes to use')
    args = parser.parse_args()
    freeze_support()
    main(args)
