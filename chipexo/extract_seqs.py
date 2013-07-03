import os,sys


def parseArgs():
    import argparse
    parser = argparse.ArgumentParser(description="Extract fasta sequences based on the names of the sequences.")
    parser.add_argument('-n', '--names', required=True, help="The file with the names of the sequences.")
    parser.add_argument('fasta', help="The fasta file")
    parser.add_argument('out', help="File to write fasta output.")
    args = parser.parse_args()
    return args


def getNames( nf ):
    f = open(nf)
    names = set()
    for r in f:
        r = r.strip()
        if r != '':
            names.add(r)
    f.close()
    return names


def main( args ):
    names = getNames( args.names )
    f = open(args.fasta)
    out = open( args.out, 'w' )
    write = False
    for r in f:
        if r[0] == '>':
            r = r.strip()
            if r[1:] in names:
                out.write(r)
                out.write('\n')
                write = True
        else:
            if write:
                out.write(r)
                write = False
    out.close()
    f.close()


if __name__=='__main__':
    args = parseArgs()
    main( args )
