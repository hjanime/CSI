import os,sys
import re
import argparse

MAX_motif = 'CAC[AG]TG[GCT]'
MAX_motif_rc = '[CGA]CA[CT]GTG'
MAX_LEN = 7
FOX_motif = '[GAT]TAAA[CT]A'
FOX_motif_rc = 'T[GA]TTTA[CTA]'
FOX_LEN = 7



def getFastaRecord( reader ):
    header = reader.readline()
    if header == '':
        return '', -1, -1, ''
    header = header[1:]
    chrom, start, end = header.strip().split('_') #1-based coordinates, inclusive
    start = int(start)
    end = int(end)
    seq = reader.readline().strip().lower()
    return chrom, start, end, seq

def parseArgs():
    parser = argparse.ArgumentParser("Get MAX and FOX motifs for the fasta records.")
    parser.add_argument('fastaFiles', nargs = '+', help="The fasta files to process")

    args = parser.parse_args()

    return args

def outputMotifCenter(chrom, start, end, seq, indexes, motif_len, withOther, out):
    '''
    Write the output for motif centers.
    '''
    only = True
    if len( indexes ) > 1:
        only  = False
    for i in indexes:
        currStart = start + abs(i) + motif_len / 2
        out.write('%s\t%d\t%d\t%s\t%s\t%s\n' % 
                ( chrom, currStart, currStart, seq[abs(i):abs(i)+motif_len], only, withOther) ) 

def makeUpper( start, end ,seq ):
    '''
    Change the specified region of the sequence to upper case.
    start: 0-based index
    end: 1-based index (+1 of the actual end)
    '''
    return seq[:start] + seq[start:end].upper() + seq[end:]


def createOutput( fastaFile ):
    '''
    Create the output files for the summary file, the motif-centered
    files for MAX and FOX motifs.
    The motif-centered files will in the following format:
    chrom <tab> start <tab> end <tab> motif <tab> onlyInSeq <tab> withOther

    onlyInSeq: indicates whether it is the only location with the motif in
               the sequence.
    withOther: indicates whether there is any locations with the other motif
               in the same sequence.
    '''
    tokens = fastaFile.split('.')
    tokens[-1] = 'txt'
    summary = open('.'.join(tokens), 'w')
    tokens[-1] = 'max_motif.txt'
    max_out = open('.'.join(tokens), 'w')
    tokens[-1] = 'fox_motif.txt'
    fox_out = open('.'.join(tokens), 'w')
    summary.write('#Chrom\tstart\tend\tsequence\tMAX_pos\tFox_pos\n')
    max_out.write('#chrom\tstart\tend\tmotif\tonlyInSeq\twithOther\n')
    fox_out.write('#chrom\tstart\tend\tmotif\tonlyInSeq\twithOther\n')
    return summary, max_out, fox_out

def writeData( chrom, start, end, seq, max_indexes,  fox_indexes, summary, max_out, fox_out):
    maxPoses = ','.join([str(i) for i in max_indexes])
    if maxPoses == '':
        maxPoses = '.'
    foxPoses = ','.join([str(i) for i in fox_indexes ])
    if foxPoses == '':
        foxPoses = '.'
    summary.write('%s\t%d\t%d\t%s\t%s\t%s\n'%(chrom, start, end, seq, maxPoses , foxPoses ))

    max_withOther = False
    fox_withOther = False
    if len(fox_indexes) > 0:
        max_withOther = True
    if len(max_indexes) > 0:
        fox_withOther = True

    outputMotifCenter( chrom, start, end, seq, max_indexes, MAX_LEN, max_withOther, max_out )
    outputMotifCenter( chrom, start, end, seq, fox_indexes, FOX_LEN, fox_withOther, fox_out )


def process( fastaFile ):
    f = open( fastaFile )
    tokens = fastaFile.split('.')
    tokens[-1] = 'txt'
    summary, max_out, fox_out = createOutput( fastaFile )
    chrom, start, end, seq = getFastaRecord( f )
    maxCount = 0
    maxSeqs = 0
    foxCount = 0
    foxSeqs = 0
    while chrom != '':
        max_indexes = []
        max_motifs = []
        max_motifs_rc = []
        fox_indexes = []
        fox_motifs = []
        fox_motifs_rc = []
        max_indexes_rc = []
        fox_indexes_rc = []
        max_search_iter = re.finditer( MAX_motif, seq, re.I )
        fox_search_iter = re.finditer( FOX_motif, seq, re.I )
        for i in max_search_iter:
            seq = makeUpper( i.start(), i.end(), seq )
            max_indexes.append( i.start() )
            max_motifs.append( seq[i.start():i.end()] )
        for i in fox_search_iter:
            seq = makeUpper( i.start(), i.end(), seq )
            fox_indexes.append( i.start() )
            fox_motifs.append( seq[i.start():i.end()] )

        max_search_iter_rc = re.finditer( MAX_motif_rc, seq, re.I )
        fox_search_iter_rc = re.finditer( FOX_motif_rc, seq, re.I )
        for i in max_search_iter_rc:
            if i.start() + 1 not in max_indexes and i.start() - 1 not in max_indexes:
                seq = makeUpper( i.start(), i.end(), seq )
                max_indexes_rc.append( -i.start() )
                max_motifs_rc.append( seq[ i.start(): i.end()] )
        for i in fox_search_iter_rc:
            seq = makeUpper( i.start(), i.end(), seq )
            fox_indexes_rc.append( -i.start() )
            fox_motifs_rc.append( seq[i.start():i.end()] )
        if len(max_indexes) + len(max_indexes_rc) > 0:
            maxSeqs += 1
        if len(fox_indexes) + len(fox_indexes_rc) > 0:
            foxSeqs += 1
        if len(max_indexes) + len(max_indexes_rc) + len(fox_indexes) + len(fox_indexes_rc) > 0:
            writeData( chrom, start, end, seq, max_indexes + max_indexes_rc,
                fox_indexes + fox_indexes_rc, summary, max_out, fox_out )
        chrom, start, end, seq = getFastaRecord( f )
    print maxSeqs,' ', foxSeqs
    summary.close()
    max_out.close()
    fox_out.close()

def main():
    args = parseArgs()
    for fasta in args.fastaFiles:
        process( fasta )



if __name__=="__main__":
    main()
