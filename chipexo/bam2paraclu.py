import pysam
import argparse
import sys

def writeToFile(items, out, chrom, strand):
    '''
    Parameters: a dictionary of the position-value tuple and a output stream.
    and the chromosome and strand of the data.
    '''
    items.sort(key = lambda k: k[0])
    for it in items:
        out.write("%s\t%s\t%d\t%d\n"%(chrom, strand, it[0], it[1]))

def main(args):
    useRegion = False
    rchrom = None
    rstart = None
    rend = None
    if args.region:
        useRegion = True
        rchrom, rstart, rend = args.region.strip().split(',')
        rstart = int(rstart)
        rend = int(rend)
    samfiles = []
    references = []
    for bamfile in args.bamfiles:
        rmode = 'r'
        if bamfile.endswith('bam'):
            rmode = 'rb'
        elif bamfile.endswith('sam'):
            rmode = 'r'
        else:
            sys.exit()
        samfiles.append(pysam.Samfile(bamfile, rmode))
        references += samfiles[-1].references
    references = set(references)
    if useRegion:
        references = [rchrom,]
    outPrefix = args.prefix
    if useRegion:
        outPrefix += '_' +args.region.replace(',','_')
    outF = open(outPrefix + "_+.para",'w')
    outR = open(outPrefix + "_-.para", 'w')
    for chrom in references:
        chromF = {}
        chromR = {}
        for samfile in samfiles:
            if not useRegion:
                chromReads = samfile.fetch(chrom)
            else:
                chromReads = samfile.fetch(chrom, rstart, rend)
            for r in chromReads:
                if r.is_reverse:
                    pos = r.positions[-1]
                    if pos in chromR:
                        chromR[pos] += 1
                    else:
                        chromR[pos] = 1
                else:
                    pos = r.pos
                    if pos in chromF:
                        chromF[pos] += 1
                    else:
                        chromF[pos] = 1
        writeToFile(chromF.items(), outF, chrom, '+')
        writeToFile(chromR.items(), outR, chrom, '-')

    for samfile in samfiles:
        samfile.close()
    outF.close()
    outR.close()


if __name__=='__main__':
    parser = argparse.ArgumentParser("Convert Bam file to input file for paraclu.")
    parser.add_argument('bamfiles',nargs='+', help="The bam file to be processed.")
    parser.add_argument('-p','--prefix', required=True, help="The prefix for the output.")
    parser.add_argument('-r', '--region', help="Only the region will be output. The region should be in format of chrom,start,end")
    args = parser.parse_args()
    main(args)
