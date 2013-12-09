import os,sys
import numpy as np


def getArgs():
    import argparse
    parser = argparse.ArgumentParser("Add motif information to the file that the fasta sequence is generated from.")
    parser.add_argument("fimoOut", help="The fimo.txt file in the fimo output folder.")
    parser.add_argument("originFile", help="The file that the motif information to add to.")
    parser.add_argument("-o", "--output", help="The file to write the output to.")

    args = parser.parse_args()
    return args

class FIMORecord:
    '''
    The FIMO record.
    Fields:
        1. Pattern_name - The name of the pattern searched.
        2. sequence_name - The name of the sequence searched against.
        3. start - relative starting position. (0-based)
        4. stop - relative end position. (1-based)
        5. strand - The strand of the hit. ('+' or '-')
        6. score - The score of the hit.
        7. pvalue - p-value
        8. qvalue - q-value
        9. matchedSeq - The matched sequence.
    '''
    def __init__(self, pattern_name, sequence_name, start, stop, strand, score, pvalue, qvalue, matchedSeq):
        self.pattern_name = pattern_name
        self.sequence_name = sequence_name
        self.start = int(start) - 1 #convert the fimo.txt 1-based coordinate into 0-based.
        self.stop = int(stop)
        self.strand = strand
        self.score = float(score)
        self.pvalue = float(pvalue)
        self.qvalue = float(qvalue)
        self.matchedSeq = matchedSeq

    def getChromLocation(self):
        '''
        Get the location information from the sequence_name.
        The sequence name is in the format:
            chrom_start_end_property,
        where start is 0-based ane end is 1-based.
        Return:
            (chrom, start, end)
        '''
        #print self.sequence_name
        chrom,start,end = self.sequence_name.split('_')[:3]
        start = int(start)
        end = int(end)
        return (chrom,start,end)

    def getBetter(self, other):
        '''
        This function is to assess whether two motifs are overlapping and if so which one to choose.
        '''
        if self.sequence_name != other.sequence_name or self.pattern_name != other.pattern_name or min(self.stop, other.stop) <= max(self.stop, other.stop):
            return [self, other]
        elif self.pvalue < other.pvalue:
            return self

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.pattern_name == other.pattern_name and \
                    self.sequence_name == other.sequence_name and \
                    self.start == other.start and \
                    self.stop == other.stop and \
                    self.score == other.score and \
                    self.pvalue == other.pvalue and \
                    self.qvalue == other.qvalue and \
                    self.matchedSeq == other.matchedSeq
        return False

    def __ne__(self, other):
        return not self.__eq__(other)


def makeMotifString(motifs):
    strings = []
    if len(motifs) > 0:
        motifs.sort(key=lambda k:(k.pattern_name, k.start, k.stop))
        curr_pattern = motifs[0].pattern_name
        curr_width = motifs[0].stop - motifs[0].start
        curr_locations = []
        curr_pvalues = []

        def createString(locs, pvals):
            locs = [str(i) for i in locs]
            pvals = ["%.2f"%(-np.log10(p)) for p in pvals]
            return "\t".join([curr_pattern, str(curr_width), str(len(locs)), ';'.join(locs), ';'.join(pvals)])

        for m in motifs:
            if m.pattern_name == curr_pattern:
                curr_locations.append(m.start)
                curr_pvalues.append(m.pvalue)
            else:
                strings.append(createString(curr_locations, curr_pvalues))
                curr_pattern = m.pattern_name
                curr_width = m.stop - m.start
                curr_locations = [m.start]
                curr_pvalues = [m.pvalue]
        strings.append(createString(curr_locations, curr_pvalues))
    return '\t'.join(strings)

def generateLargeTableHeader(numMotifs):
    tokens = ["#chrom",
              "start",
              "end",
              "name",
              "score",
              "strand",
              "signalValue",
              "pValue",
              "qValue",
              "summit",
              "shift",
              "super-enhancer_id",
              "type",
              "genes"]
    motif_headers = ["pattern", "pattern_width", "#hits", "hitStarts", "hitPvalues"]
    for i in range(numMotifs):
        tokens += motif_headers
    return '\t'.join(tokens)

def loadOrigin(originFile):
    '''
    Load data from the file to add the motifs to.
    Assume that the first three columns of the origin file contains chrom, start and end.
    And the format of the coordinates is assumed as the same as BED format.
    Parameters:
        originFile - the filename of the origin file.
    Output:
        {chrom:{start:end, start1:end1...}, chrom1:{start2:end2, start3:end3 ...}}
    '''
    origin = {}
    with open(originFile) as f:
        for r in f:
            if r.startswith('#'):
                continue
            chrom, start, end = r.strip().split('\t')[:3]
            start = int(start)
            end = int(end)
            if chrom not in origin:
                origin[chrom] = {}
            origin[chrom][start] = [end, r.strip()]
    return origin

def loadFIMOOut(fimoOut):
    '''
    Load the fimo output.
    return
    '''
    motifs = []
    patterns = set()
    with open(fimoOut) as f:
        for r in f:
            if r.startswith('#'):
                continue
            motifs.append(FIMORecord(*r.strip().split('\t')))
            if motifs[-1].pattern_name not in patterns:
                patterns.add(motifs[-1].pattern_name)
    motifs.sort(key=lambda k:(k.sequence_name, k.pattern_name, k.start, k.stop))
    motifs_nodup = []
    if len(motifs) > 0:
        prev = motifs[0]
        motifs_nodup.append(prev)
        for m in motifs[1:]:
            if m == prev:
                continue
            motifs_nodup.append(m)
            prev = m
    return motifs_nodup, patterns

def addMotifToOrigin(motifs, origin):
    for m in motifs:
        chrom,start,end = m.getChromLocation()
        if chrom in origin and start in origin[chrom] and origin[chrom][start][0] == end:
            origin[chrom][start].append(m)

def writeData(out, origin):
    for chrom in origin:
        for start in origin[chrom]:
            curr = origin[chrom][start][1:] #the first element is the end, omit it.
            motif_string = makeMotifString(curr[1:])
            out.write('%s\t%s\n'%(curr[0], motif_string))




def main(fimoOut, originFile, output=None):
    fimo, patterns = loadFIMOOut(fimoOut)
    origin = loadOrigin(originFile)
    addMotifToOrigin(fimo, origin)
    if output is None:
        output = originFile + ".withMotif.xls"
    with open(output, 'w') as out:
        out.write(generateLargeTableHeader(len(patterns)))
        out.write('\n')
        writeData(out, origin)


if __name__=='__main__':
    args = getArgs()
    main(args.fimoOut, args.originFile, output=args.output)
