import os,sys


def loadData(fh):
    seqs = {}
    count = 0
    for r in fh:
        if r[0] == '#':
            continue
        tokens = r.strip().split()
        if tokens[1] not in seqs:
            seqs[tokens[1]] = set()
        tempStrand = '.'
        if tokens[4] == '+':
            tempStrand = '-'
        elif tokens[4] == '-':
            tempStrand = '+'
        temp = (tokens[0],tokens[2],tokens[3],tempStrand,tokens[-1])
        if temp not in seqs[tokens[1]]:
            toAdd = (tokens[0],tokens[2],tokens[3],tokens[4],tokens[-1])
            seqs[tokens[1]].add(toAdd)
            count += 1
    return seqs, count

def getSeqsWithOneHit(seqs):
    one_hit_seqs = []
    for k in seqs:
        if ken(seqs

def main(inputfile):
    with open(inputfile) as fh:
        seqs, total_hits = loadData(fh)
        total_seqs = len(seqs.keys())
        seqs_with_only_one_hit = 0






def getArgs():
    import argparse

    parser = argparse.ArgumentParser(description="Get the motif statistics from the fimo txt output.")
    parser.add_argument("inputfile", help="The input file to process. fimo.txt output of fimo.")
    args = parser.parse_args()
    return args

if __name__=='__main__':
    args = getArgs()
    main(*args)
