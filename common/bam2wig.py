import argparse
import os,sys
from contextlib import contextmanager, closing
import logging
from collections import defaultdict
import re


import pysam
import numpy as np
logging.basicConfig(level=logging.INFO)
info = logging.info

def main(bamfile, prefix=None, chrom='all',
         start=0, end=None, thresh=0):
    if not os.path.exists(bamfile):
        sys.stderr.write("Bam file not found. Exiting.")
        sys.exit(1)
    if prefix is None:
        prefix = os.path.splitext(bamfile)[0]
    if start != 0:
        start = int(start) - 1
    if end is not None:
        end = int(end)

    regions = [(chrom,start,end)]
    wig_f = '%s_Forward.wig'%prefix
    wig_r = '%s_Reverse.wig'%prefix
    with indexed_bam(bamfile) as reader:
        with open(wig_f,'w') as outF:
            with open(wig_r,'w') as outR:
                sizes = zip(reader.references, reader.lengths)
                if len(regions) == 1 and regions[0][0] == 'all':
                    regions = [(name,0,length) for name,length in sizes if re.match('^chr(\d+|X|Y|M)$', name)]
                writeTrackLine(outF, bamfile)
                writeTrackLine(outR, bamfile)
                loadData(regions, reader, outF, outR, thresh)

@contextmanager
def indexed_bam(bamfile):
    if not os.path.exists(bamfile + '.bai'):
        pysam.index(bamfile)
    reader = pysam.Samfile(bamfile,'rb')
    yield reader
    reader.close()

def writeTrackLine(out_handle, bamfile):
    out_handle.write("track %s\n" % " ".join(["type=wiggle_0",
        "name=%s" % os.path.splitext(os.path.split(bamfile)[-1])[0],
        "visibility=full",
        ]))

def loadData(regions, reader, outF, outR, thresh):
    for (chrom,start,end) in regions:
        info("Processing %s---"%chrom)
        if end is None and chrom in reader.references:
            end = reader.lengths[reader.references.index(chrom)]
        tempFcount = defaultdict(int)
        tempRcount = defaultdict(int)
        count = 0
        for r in reader.fetch(chrom, start, end):
            if not r.is_reverse:
                tempFcount[r.pos] += 1
            else:
                if r.aend is not None:
                    tempRcount[r.aend - 1] += 1
                elif r.pos + r.qlen - 1 < end:
                    tempRcount[r.pos + r.qlen - 1] += 1
                else:
                    tempRcount[end - 1] += 1
            count += 1
        info("Total %d"%count)
        writeData(tempFcount, outF, chrom, start, thresh)
        writeData(tempRcount, outR, chrom, start, thresh)


def writeData(data, handle, chrom, start, thresh):
    handle.write('variableStep chrom=%s\n'%chrom)
    outData = []
    keys = data.keys()
    keys.sort()
    for k in keys:
        d = data[k]
        if d > thresh:
            outData.append('%d\t%.1f\n'%(start+k+1, d))
            #handle.write('%s\t%.1f\n'%(start+i+1, d))
    handle.write(''.join(outData))



if __name__=='__main__':
    parser = argparse.ArgumentParser("Convert bam to wiggle files separated by strand.")
    parser.add_argument('bam', help="The bam input.")
    parser.add_argument('-o', '--prefix', help="The output prefix. Default: the input name without extension.")
    parser.add_argument('-c','--chrom', default="all", help="The chromosome to use. Default: use all.")
    parser.add_argument('-s','--start', type=int, help="The start position of the region of interest. 1-based.")
    parser.add_argument('-e','--end', type=int, help="The end position of the region of interest.")
    parser.add_argument('-t','--thresh', type=float, help="Below (inclusive) which the data point should be ignored.")

    args = parser.parse_args()
    kwargs = dict(prefix = args.prefix,
                  chrom = args.chrom or "all",
                  start = args.start or 0,
                  end = args.end,
                  thresh = args.thresh or 0,
                  )
    main(args.bam, **kwargs)


