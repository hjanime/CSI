import argparse

def main(args):
    f = open(args.paraout)
    out = open(args.bed, 'w')
    for r in f:
        if r[0] == '#':
            continue
        chrom,strand,start,end,numposes,dsum,mindense,maxdense = r.strip().split()
        score = dsum
        if args.score == 'numposes':
            score = numposes
        elif args.score == 'mindense':
            score = mindense
        elif args.score == 'maxdense':
            score = maxdense
        name = '_'.join(r.strip().split())
        out.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(chrom,start,end,name,score,strand))
    out.close()

if __name__=='__main__':
    parser = argparse.ArgumentParser("Convert paraclu output to bed format.")
    parser.add_argument("paraout", help="The paraclu output file.")
    parser.add_argument("bed", help="The name of the bed file to be written to.")
    parser.add_argument("-s", '--score',default='sum', choices=['numposes','sum','mindense','maxdense'],
                        help="The column to use as the bed file score.")
    args = parser.parse_args()
    main(args)
