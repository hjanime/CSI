import os,sys


def getArgs():
    import argparse
    parser = argparse.ArgumentParser(description="Add  cw_distance=xx info to the header of fasta sequence headers.")

    parser.add_argument('fasta', help="The fasta file to modify. If the input is XXX.fa, the output will be XXX.modified.fa.")
    parser.add_argument('annotation', help="The annotation file that generates the fasta file. MUST BE IN GFF or BED FORMAT.")
    parser.add_argument('distCol', type=int, help="The index of the column that contains the dist information. Index starts from 0")
    parser.add_argument('--format', choices=['BED','GFF'], default='GFF', help="The format of the annotation file. It can be either BED or GFF.")
    parser.add_argument('--header', action='store_true', default=False, help="Whether the first line should be ignored.")
    parser.add_argument('--onebased', default=False, action='store_true', help="Whether the fasta header coordinates are 1-based. Default is the same format as BED.")


    args = parser.parse_args()
    return args

def getAnnotation( annFile, distCol, hasHeader, fastaIsOneBased, annFormat ):
    f = open( annFile )
    if hasHeader:
        f.readline()
    lookupDict = {}
    for r in f:
        if r.startswith('#'):
            continue
        seqKey = None
        chrom = None
        start = None
        end = None
        tokens = r.strip().split('\t')
        if annFormat == 'GFF':
            chrom = tokens[0]
            start = int( tokens[3] )
            end = int( tokens[4] )
            if not fastaIsOneBased:
                start -= 1
        else:
            chrom = tokens[0]
            start = int( tokens[1] )
            end = int( tokens[2] )
            if fastaIsOneBased:
                start += 1
        seqKey = '%s_%d_%d'%(chrom, start, end)
        lookupDict[ seqKey ] = tokens[distCol]
    return lookupDict

def main( args ):
    import re
    ann = getAnnotation( args.annotation, args.distCol, args.header, args.onebased, args.format)
    f = open(args.fasta)
    ftokens = args.fasta.split('.')
    ftokens[-1] = 'modified.fa'
    out = open( '.'.join( ftokens ), 'w' )
    for r in f:
        if r.startswith('>'):
            header = r.strip()[1:]
            header = header.replace(':', '_')
            header = header.replace('-','_')
            if not re.match('[a-zA-Z0-9]+_\d+_\d+', header):
                print 'Header in wrong format: ' + r
                continue
            elif header not in ann:
                print 'Header not found in annotation: ' + r
                continue
            header = header + '_' + ann[ header ]
            out.write('>%s\n' % header )
        else:
            out.write( r )
    f.close()
    out.close()


if __name__=='__main__':
    args = getArgs()
    main(args)
