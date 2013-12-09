import os,sys
import argparse

def main(args):
    for filename in args.inputFiles:
        f = open(filename)
        tokens = filename.split('.')
        tokens[-1] = str(args.extend) + "_center.bed"
        out = open('.'.join(tokens),'w')
        written = {}
        for r in f:
            if r.startswith('#'):
                continue
            tokens = r.strip().split('\t')
            chrom = tokens[0]
            start = 0
            end = 0
            if len(tokens[0].split('_')) > 1:
                chrom, start, end = tokens[0].split('_')
            if chrom not in written:
                written[chrom] = set()
            start = int(start)
            end = int(end)
            currPos = ( int( tokens[3] ) + int( tokens[4] ) ) / 2
            center = start + currPos - 1
            if center not in written[chrom]:
                written[chrom].add(center)
                cstart = center - args.extend
                cend = center + args.extend
                out.write('\t'.join([chrom, str(cstart), str(cend),'.']))
                out.write('\n')
        out.close()


def parseArgs():
    parser = argparse.ArgumentParser(description = "Get the center locations from the gff files of fimo output.")
    parser.add_argument("inputFiles", nargs='+', help="The input files.")
    parser.add_argument('-e','--extend', type=int, default = 0, help="How much to extend on both side of the center.")
    args = parser.parse_args()
    return args



if __name__=='__main__':
    args = parseArgs()
    main( args )
