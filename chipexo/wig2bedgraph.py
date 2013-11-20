import wig
import argparse

def getArgs():
    parser = argparse.ArgumentParser("Convert wig file to bedgraph file.")
    parser.add_argument('filename', help="Input wig file")
    parser.add_argument('strand', choices=['+','-','.'], help="The strand of the input. Can be [+, - , . ]")
    parser.add_argument('output', help="The output file.")

    args = parser.parse_args()
    return args


if __name__=='__main__':
    args = getArgs()
    wigData = wig.loadWig(args.filename, False, args.strand)
    wig.writeAsBedGraph(wigData, args.output)
