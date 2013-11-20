def getData(array, chrom, start, end, fh):
    for l in fh:
        if l.startswith('chr8'):
            tokens = l.strip().split()
            tstart = int(tokens[1])
            tend = int(tokens[2])
            tv = float(tokens[3])
            if tstart >= start and tend <= end:
                for i in range(tstart, tend):
                    array[tstart-start] = tv

def normalize(inA, outA ):
    aveA = inA.mean()
    stdA = inA.std()
    for i in range(inA.shape[0]):
        outA[i] = (inA[i] - aveA) / stdA
