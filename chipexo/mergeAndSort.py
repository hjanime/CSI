import os,sys
import argparse

SID = { 'length': 3,
        'strand': 4,
        'summit': 5,
        'height': 6,
        'height5': 8,
        'pvalue': 9,
        'fold'  : 10,
        'qvalue': 11 }

def getRecords( filename ):
    f = open( filename )
    records = {} # store the records by chromosome to reduce the sorting time.
    for r in f:
        if r.strip() == '' or r[0] == '#':
            pass
        else:
            tokens = r.strip().split('\t')
            tokens[1] = int( tokens[1] )
            tokens[2] = int( tokens[2] )
            tokens[3] = int( tokens[3] )
            tokens[5] = int( tokens[5].split(',')[0] )
            tokens[6] = float( tokens[6] )
            tokens[7] = int( tokens[7] )
            tokens[8] = float( tokens[8])
            tokens[9] = float( tokens[9] )
            tokens[10] = float( tokens[10] )
            tokens[11] = float( tokens[11] )
            tokens[-1] = float( tokens[-1] )
            tokens[-2] = int(tokens[-1])
            if tokens[0] not in records:
                records[ tokens[0] ] = []
            records[ tokens[0] ].append( tokens )
    f.close()
    return records

def sortRecords( records, sortIdx ):
    for key in records:
        '''
        if sortType == 'position': #sort by peak location
            records[ key ].sort( key = lambda k: ( int( k[1] ), int( k[2] ) ) )
        elif sortType == 'length': #sort by peak length
            records[ key ].sort( key = lambda k: ( int( k[3] ) ) )
        elif sortType == 'summit': #sort by summit location
            records[ key ].sort( key = lambda k: ( int( k[5] ) ) )
        elif sortType == 'height': #sort by summit height
            records[ key ].sort( key = lambda k: ( float( k[6] ) ) )
        elif sortType == 'fold': #sort by fold enrichment
            records[ key ].sort( key = lambda k: ( float( k[8] ) ) )
        elif sortType == 'pvalue': #sort by log10(pvalue)
            records[ key ].sort( key = lambda k: (float( k[7] ) ) )
        elif sortType == 'qvalue': #sort by log10(qvalue)
            records[ key ].sort( key = lambda k: (float( k[9] ) ) )
        '''
        records[ key ].sort( key = lambda k: tuple( [ k[i] for i in sortIdx ] ) )


def writeRecords( records, outfile, header ):
    out = open( outfile, 'w' )
    if header != None:
        out.write( header )
    chroms = records.keys()
    chroms.sort()
    for chrom in chroms:
        for r in records[chrom]:
            out.write( '\t'.join( [ str(i) for i in r ] ) )
            out.write( '\n' )
    out.close()

def getHeader( filename ):
    f = open( filename )
    header = []
    for r in f:
        if r[0] == '#':
            header.append( r )
        else:
            break
    f.close()
    return ''.join( header )

if __name__ == '__main__':
    parser = argparse.ArgumentParser( description = "Merge files and sort the records. Note: the commented lines will all be discarded. Comment lines start with '#'." + 
                                                    "The files used by this program must be the xls output of getSummit.py. For other file fomat such as BED or GFF, " +
                                                    "Bedtools should be used for convenience.")
    parser.add_argument( '-i', '--inputs', required = True, nargs = '+', help = "The input names." )
    parser.add_argument( '-o', '--output', required = True, help="The output file name." )
    parser.add_argument( '-s', '--sortType', nargs = '*', default = ['position',], help = "The keys used for sorting."
                + " \nMultiple keys could be used and the"
                + " \nfile will be sorted by the order of the keys."
                + " \nIn any case, the records will always be sorted by chromosome first."
                + " \nThe values can be: position, length, summit," 
                + " \nheight, strand, pvalue, fold, qvalue")
    parser.add_argument( '-H', '--header', action = "store_true", default = False, help = "Whether to print the header of the first file in the list to the output." )

    args = parser.parse_args()

    print args

    header = None
    if args.header:
        header = getHeader( args.inputs[0] )

    records = {}
    for i in args.inputs:
        tempRecords = getRecords( i )
        for key in tempRecords:
            if key not in records:
                records[ key ] = tempRecords[ key ]
            else:
                records[ key ] += tempRecords[ key ]
    sortIdx = []
    if args.sortType:
        for s in args.sortType:
            if s.lower() == "position":
                sortIdx += [1,2,]
            else:
                sortIdx.append( SID[s.lower()] )
    sortRecords( records, sortIdx )

    writeRecords( records, args.output, header )

