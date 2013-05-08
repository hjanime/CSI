import os,sys
import wig
import gc

from multiprocessing import Process, Queue, current_process, freeze_support


offset = 0

def add_chrom_data(taskQ, outQ, processID ):
    for chrom, chromWig in iter( taskQ.get, 'STOP' ):
        print "add chrom data Process ", processID," is processing ", chrom
        lines = []
        lines.append("variableStep chrom=%s\n"%(chrom,))
        startp = chromWig[0,0]
        expanded = wig.expandWig( chromWig, offset, 1 )
        for i in range( expanded.shape[0] ):
            if expanded[i] > 0.1:
                lines.append( "%d\t%f\n"%(int(i + startp - offset), expanded[i], ) )
        #for i in range( chromWig.shape[0] ):
        #    lines.append('%d\t%f\n'%( chromWig[ i, 0], chromWig[i,1]))
        expanded = None
        outQ.put(lines)
        gc.collect()

def main():
    for filename in sys.argv[1:]:
        lines = []
        wigdata = wig.loadWig( filename, False )
        print wigdata
        tokens = filename.split('.')
        tokens[-1] = "smoothed.wig"

        chroms = Queue()
        count = 0
        outQ = Queue()
        for chrom in wigdata:
            print "add ", chrom
            count += 1
            chroms.put( (chrom, wigdata[chrom][:,0:2]) )
        wigdata = None

        NUM_PROCESSES = 3

        
        processID = 1
        for i in range( NUM_PROCESSES):
            Process( target=add_chrom_data, args=( chroms, outQ, processID ) ).start()
            processID += 1

        out = open( '.'.join(tokens), 'w' )

        for i in range( count ):
            out.write( ''.join(outQ.get()) )


        out.close()
        for i in range( NUM_PROCESSES ):
            chroms.put('STOP')





if __name__=='__main__':
    freeze_support()
    main()
