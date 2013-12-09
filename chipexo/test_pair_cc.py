import os,sys
import numpy as np
from operator import itemgetter
import bisect

sys.path.insert(0,'/home/caofan/Documents/scripts/rnaseq')
from readBed import BEDReader
from getSummit import loadRejectRegion, isReject
from sortedcollection import SortedCollection

def normalize(inA, outA ):
    aveA = inA.mean()
    stdA = inA.std()
    for i in range(inA.shape[0]):
        outA[i] = (inA[i] - aveA) / stdA

class Peak:
    '''
    The class should contain all the fields from narrowPeak format.

    NarrowPeak format fields:

    1. chrom - Name of the chromosome (or contig, scaffold, etc.).
    2. chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
    3. chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
    4. name - Name given to a region (preferably unique). Use '.' if no name is assigned.
    5. score - Indicates how dark the peak will be displayed in the browser (0-1000). If all scores were '0' when the data were submitted to the DCC, the DCC assigned scores 1-1000 based on signal value. Ideally the average signalValue per base spread is between 100-1000.
    6. strand - +/- to denote strand or orientation (whenever applicable). Use '.' if no orientation is assigned.
    7. signalValue - Measurement of overall (usually, average) enrichment for the region.
    8. pValue - Measurement of statistical significance (-log10). Use -1 if no pValue is assigned.
    9. qValue - Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned.
    10. peak - Point-source called for this peak; 0-based offset from chromStart. Use -1 if no point-source called.

    In addition:
    11. data - an 2-d array of the values at each location within the peak. data[0] is the forward strand and data[1] is the reverse strand.
    12. ctrl - similar to data field, it stores the control data at the peak.
    13. processed - a bool field to track the status of the current peak. It can be used for any purposes.

    The class needs a tuple/list ordered as described to instantiate.
    The tuple/list should have at least fields 1-3.
    '''
    def __init__(self, chrom, start, end, name = '.', score = 0, strand = '.', signalValue = None, pValue = -1, qValue = -1, summit = -1):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = name
        self.score = score
        self.strand = strand
        self.signalValue = signalValue
        self.pValue = pValue
        self.qValue = qValue
        self.summit = summit

        self.data = None
        self.ctrl = None
        self.cc = None
        self.shift = None
        self.old_shift= None
        self.processed = False
        self.size = self.end - self.start
        self.stitchedId = 'NA'
        if not (isinstance(self.start, int) and isinstance(self.end, int)):
            raise ValueError("start and end should be integers.\n")



    def setStitchID(self, idx):
        self.stitchedId = idx

    def prepareData(self, isControl):
        '''
        Prepare a numpy 1-d array for loading coverage data.
        Output:
            None. The data field will be a 1-d array of 0s with size end-start.
        '''
        if not isControl and self.data == None:
            self.data = [np.zeros(self.size, dtype='float32'), np.zeros(2*self.size-1, dtype='float32')]
        elif isControl and self.ctrl == None:
            self.ctrl = np.zeros((2, self.size), dtype='float32')


    def _addData(self, data, dstart, dend, v, strand=0):
        '''
        Helper function to perform the data adding.
        '''
        dstart = int(dstart)
        dend = int(dend)
        assert dstart < dend
        assert isinstance(dstart, int) and isinstance(dend, int) and isinstance(v, float)
        assert strand == 0 or strand == 1, "Strand is not correct"
        if self.start <= dstart < self.end or self.start < dend <= self.end:
            for i in range(max(dstart,self.start), min(self.end, dend)):
                data[strand][i-self.start] += v
        elif strand == 1 and self.start <= dstart < self.end + self.size - 1 or self.start < dend <= self.end + self.size - 1:
            for i in range(max(dstart,self.start), min(self.end+self.size-1, dend)):
                data[strand][i-self.start] += v

    def addData(self, dstart, dend, v ,strand = 0):
        '''
        Add value to the data array.
        Parameters:
            dstart: the start of the data. Absolute coordinate on the chromosome
                    Integer, 0-based.
            dend: the end of the data. Absolute coordinate on the chromosome.
                  Integer, 1-based, not inclusive.
            v: the value of the segment to be added. Float.
            strand: optional, default to 0. The index of the data array to add to. Can only be 0 or 1.
        '''
        self._addData(self.data, dstart, dend, v, strand)
        self.processed = False #data modified

    def addCtrl(self, dstart, dend, v, strand = 0):
        '''
        Add value to the ctrl array.
        Parameters:
            dstart: the start of the data. Absolute coordinate on the chromosome
                    Integer, 0-based.
            dend: the end of the data. Absolute coordinate on the chromosome.
                  Integer, 1-based, not inclusive.
            v: the value of the segment to be added. Float.
            strand: optional, default to 0. The index of the data array to add to. Can only be 0 or 1.
        '''
        self._addData(self.ctrl, dstart, dend, v, strand)

    def getCC(self, smooth=False):
        '''
        Calculate the cross-correlation of two peaks.
        Output:
            cc: the cross-correlation using np.correlate(peak.data, self.data, 'full').
                Like setting the other peak stable and move the current peak around.
        '''
        assert self.data != None, "Data is empty"
        if self.cc != None and self.processed:
            return self.cc

        '''
        data1 = self.data[0]
        data2 = self.data[1]
        data2 = (self.data[1] - self.data[1].mean())/self.data[1].std()
        data1 = (self.data[0] - self.data[0].mean())/self.data[0].std()
        if smooth:
            w = np.blackman(11)
            data1 = np.convolve(w, self.data[0],mode='same')
            data2 = np.convolve(w, self.data[1],mode='same')
        #size1 = peak.end - peak.start
        #size2 = self.end - self.start
        #ratio = max(size1, size2) * 1.0 / min(size1, size2)
        #ccn = ccn/ratio
        #x = np.argmax(ccn)
        #shift = peak.start + x + 1 - self.end

        self.cc = np.correlate(data2, data1, 'full')
        ccn = np.zeros_like(self.cc)
        size = data1.shape[0]
        for i in range(len(self.cc)):
            ccn[i] = self.cc[i] / (min(i + 1, size) - max(i + 1 - size, 0))
        for i in range(len(self.cc)):
            if i < data1.shape[0]:
                self.cc[i] = self.cc[i] / (i+1)
            else:
                self.cc[i] = self.cc[i] / (2*data1.shape[0] - i - 1)
        self.shift = np.argmax(ccn[data1.shape[0]-1:data1.shape[0]*3/2])
        '''
        w = np.blackman(11)
        data1 = self.data[0]
        data2 = self.data[1]
        mean1 = data1[:self.size].mean()
        mean2 = data2[:self.size].mean()
        data1 = data1 / mean1
        data2 = data2 / mean2
        data1 = np.convolve(w, self.data[0],mode='same')
        data2 = np.convolve(w, self.data[1],mode='same')
        #data1 = data1*np.sqrt(data1)
        #data2 = data2*np.sqrt(data2)
        maxI = 0
        maxScore = 0
        cc = np.correlate(data2, data1[:self.size], 'valid')
        #self.shift= np.argmax(cc[data1.shape[0]-1:])
        self.shift= np.argmax(cc)
        individual = data1[:self.size]*data2[self.shift:self.shift+self.size]
        self.summit = self.shift / 2 + np.argmax(individual)

        self.processed = True
        self.cc = cc
        return self.cc

    def getOldShift(self):
        w = np.blackman(11)
        data1 = self.data[0][:self.size]
        data2 = self.data[1][:self.size]
        mean1 = data1.mean()
        mean2 = data2.mean()
        data1 = data1 / mean1
        data2 = data2 / mean2
        data1 = np.convolve(w, self.data[0],mode='same')
        data2 = np.convolve(w, self.data[1],mode='same')
        data1[data1<=0]=0.001
        data2[data2<=0]=0.001
        data1 = data1*np.sqrt(data1)
        data2 = data2*np.sqrt(data2)
        maxI = 0
        maxScore = 0
        cc = np.correlate(data2, data1, 'full')
        self.old_shift = np.argmax(cc[data1.shape[0]-1:])


    def __str__(self):
        return '\t'.join([str(i) for i in [self.chrom, self.start, self.end, self.name, self.score, self.strand, self.signalValue, self.pValue, self.qValue, self.summit, self.shift, self.stitchedId]])

    def simpleStr(self):
        a = '%s_%d_%d_%f'%(self.chrom, self.start, self.end, self.signalValue)
        if self.shift != None:
            a = '%s_%d'%(a, self.shift)
        return a

    def plot(self,savefig=None):
        '''
        '''
        import pylab as pl
        if self.end-self.start < 51:
            return
        x = np.arange(self.start, self.end)
        xr = np.arange(self.start, self.end - self.start + self.end - 1)
        pl.figure()
        pl.subplot(221)
        pl.plot(x, self.data[0])
        pl.plot(xr, self.data[1], color='r')
        pl.subplot(223)
        pl.plot(x+self.shift/2, self.data[0])
        pl.plot(xr-(self.shift - self.shift/2), self.data[1], color='r')
        #pl.subplot(325)
        #pl.plot(x+self.old_shift/2, self.data[0])
        #pl.plot(xr-(self.old_shift - self.old_shift/2), self.data[1], color='r')

        w = np.blackman(51)
        data1 = np.convolve(w, self.data[0],mode='same')
        data2 = np.convolve(w, self.data[1],mode='same')
        pl.subplot(222)
        pl.plot(x, data1)
        pl.plot(xr, data2, color='r')
        pl.subplot(224)
        pl.plot(x+self.shift/2, data1)
        pl.plot(xr-(self.shift - self.shift/2), data2, color='r')
        #pl.subplot(326)
        #pl.plot(x+self.old_shift/2, data1)
        #pl.plot(xr-(self.old_shift - self.old_shift/2), data2, color='r')


        if savefig != None:
            pl.savefig(savefig+"/"+self.simpleStr()+'_'+str(self.old_shift)+'.png')
            pl.cla()
            pl.clf()
            pl.close()
        else:
            pl.show()




class StitchedPeaks:
    '''
    Class for stiched peaks.
    Fields:
        1. peaks: a list of the peaks of the Peak class.
        2. stichDist: the distance within which peaks are stitched.
    '''
    def __init__(self, index, peaks, stitchDist):
        self.index = index
        self.peaks = []
        self.peaks += peaks
        for p in self.peaks:
            p.setStitchID(self.index)
        self._sort()
        self.stitchDist = stitchDist
        self.totalSignalInConst = -1
        self.totalSignal = -1
        self.totalCtrlInConst = -1
        self.totalCtrl = -1

    def getChrom(self):
        if len(self.peaks) <= 0:
            return None
        return self.peaks[0].chrom

    def getStart(self):
        if len(self.peaks) <= 0:
            return None
        return self.peaks[0].start

    def getEnd(self):
        if len(self.peaks) <= 0:
            return None
        return self.peaks[-1].end

    def _calDataInConst(self, control=False):
        total = 0
        if not control:
            for p in self.peaks:
                if p.data != None:
                    total += p.data[0].sum() + p.data[1][:p.size].sum()
        else:
            for p in self.peaks:
                if p.ctrl != None:
                    total += p.ctrl.sum()
        return total



    def getTotalSignalInConst(self):
        if self.totalSignalInConst != -1:
            return self.totalSignalInConst
        self.totalSignalInConst = self._calDataInConst()
        return self.totalSignalInConst

    def getTotalCtrlInConst(self):
        if self.totalCtrlInConst != -1:
            return self.totalCtrlInConst
        self.totalCtrlInConst = self._calDataInConst(True)
        return self.totalCtrlInConst

    def calTotalSignal(self, fwig, rwig):
        import bisect
        if fwig == None or rwig == None:
            return -1
        def _getBoundIdx(iwig, chrom, tstart, tend):
            s = bisect.bisect_left(iwig[chrom][:,0], tstart)
            e = bisect.bisect_right(iwig[chrom][:,0], tend)
            return s,e

        if len(self.peaks) <= 0:
            return -1
        total = -1
        if self.getChrom() in fwig and self.getChrom() in rwig:
            chrom = self.getChrom()
            start = self.getStart()
            end = self.getEnd()
            fs, fe = _getBoundIdx(fwig, chrom, start, end)
            rs, re = _getBoundIdx(rwig, chrom, start ,end)
            total = fwig[chrom][fs:fe,1].sum() + rwig[chrom][rs:re,1].sum()
        return total

    def getTotalSignal(self, fwig, rwig):
        if self.totalSignal != -1:
            return self.totalSignal
        self.totalSignal = self.calTotalSignal(fwig, rwig)
        return self.totalSignal

    def getTotalCtrl(self, fwig, rwig):
        if self.totalCtrl != -1:
            return self.totalCtrl
        self.totalCtrl = self.calTotalSignal(fwig, rwig)
        return self.totalCtrl


    def _sort(self):
        self.peaks.sort(key=lambda k:(k.chrom, k.start, k.end))

    def addPeak(self, p):
        '''
        Add a peak to the list.
        '''
        self.peaks.append(p)
        self._sort()

    def addPeaks(self, ps):
        '''
        Add a list of peaks.
        '''
        self.peaks += ps
        self._sort()

    def __str__(self):
        a = []
        a.append(self.getChrom())
        a.append(self.getStart())
        a.append(self.getEnd())
        a.append(self.peaks[0].name+'_stitched')
        a.append(self.index)
        a.append(len(self.peaks))
        a.append(self.totalSignal)
        a.append(self.totalCtrl)
        a.append(self.totalSignalInConst)
        a.append(self.totalCtrlInConst)
        a.append(','.join([p.simpleStr() for p in self.peaks]))
        return '\t'.join([str(i) for i in a])


def getStitchHeader():
    a = ['#Chrom',
         'Start',
         'End',
         'Name',
         'ID',
         'Num_peaks',
         'Total_signal',
         'Total_ctrl',
         'Total_sig_in_const',
         'Total_ctrl_in_const',
         'Consts']
    return '\t'.join(a)


class PeakPair:
    '''
    Class to store peak pairs.
    Fields:
        1. fpeak: a Peak object, forward peak.
        2. rpeak: a Peak object, reverse peak. During initiation, the class will determine
                  direction of the two peaks if any of them contain strand information.
                  If no strand information found or same strands found, the order of initiation
                  will be followed.
        3. ccn: The cross-correlation profile normalized by the length of the overlapping region.
        4. shift: the shift between the two peaks.
    '''
    def __init__(self, peak1, peak2):
        #assert isinstance(peak1, Peak) and isinstance(peak2, Peak), "Wrong type"
        self.fpeak = peak1
        self.rpeak = peak2
        if (self.fpeak.strand == '-' and self.rpeak.strand != '-') or (self.rpeak.strand == '+' and self.fpeak.strand != '+'):
            self.fpeak = peak2
            self.rpeak = peak1
        self.ccn = None
        self.shift = None
        if self.fpeak.data != None and self.rpeak.data != None:
            self.ccn = self.fpeak.getCC(self.rpeak)
            #for i in range(self.fpeak.end - self.rpeak.start - 1):
            #    self.ccn[i] *= -1
            self.shift = np.argmax(self.ccn) + self.rpeak.start + 1 - self.fpeak.end



def loadNarrowPeaks(reader, minsize=0):
    '''
    Parameters:
        reader: the BEDReader to read in the narrowPeak files
    Output:
        peaks: {chrom: [Peak, Peak...], chrom2: [Peak, Peak...]...}
    '''
    peaks = {}
    for r in reader:
        if r[2] - r[1] < minsize:
            continue
        if r[0] not in peaks:
            peaks[r[0]] = [Peak(*r),]
        else:
            peaks[r[0]].append(Peak(*r))
    for chrom in peaks:
        peaks[chrom].sort(key=lambda k:(k.start, k.end,))
    return peaks

def addWigToPeak(fw, rw, peaks, isControl=False):
    '''
    Parameters:
        fw: the wig file of forward strand.
        rw: the wig file of reverse strand.
        peaks: the peaks to add the data to.
    output:
        None. The peaks will be populated with data.
    '''
    import wig
    fwig = wig.loadWig(fw, False, '+')
    rwig = wig.loadWig(rw, False, '-')
    total_f = 0
    total_r = 0
    chroms = []
    for chrom in peaks:
        if not (chrom in fwig and chrom in rwig):
            continue
        chroms.append(chrom)
        chromPeaks = peaks[chrom]
        chromFwig = fwig[chrom]
        chromRwig = rwig[chrom]
        #total_f += chromFwig[:,1].sum()
        #total_r += chromRwig[:,1].sum()
        fi = 0 #index to track chromFwig
        ri = 0 #index to track chromRwig
        pi = 0 #index to track chromPeaks
        while pi < len(chromPeaks):
            p = chromPeaks[pi]
            p.prepareData(isControl)
            #add the forward data

            def _addWig(cp, cwig, i1, strand):
                while i1 < cwig.shape[0] and cwig[i1,0]-1 < cp.start:
                    i1 += 1

                while i1 < cwig.shape[0] and cp.start <= cwig[i1,0] - 1 < cp.end:
                    if isControl:
                        cp.addCtrl(cwig[i1,0]-1, cwig[i1,0], cwig[i1,1], strand)
                    else:
                        cp.addData(cwig[i1,0]-1, cwig[i1,0], cwig[i1,1], strand)
                    i1 += 1

                if strand == 1 and not isControl:
                    i2 = i1
                    cend = cp.end - cp.start + cp.end - 1
                    if pi < len(chromPeaks)-1 and chromPeaks[pi+1].start < cend:
                        cend = chromPeaks[pi+1].start
                    while i2 < cwig.shape[0] and cwig[i2,0] - 1 < cend:
                        cp.addData(cwig[i2,0]-1, cwig[i2,0], cwig[i2,1], strand)
                        i2 += 1

                return i1

            fi = _addWig(p, chromFwig, fi, 0)
            ri = _addWig(p, chromRwig, ri, 1)

            pi += 1
    return fwig, rwig, chroms

def filterPeaks(peaks, num_poses, rpm, ratio, chroms, rejected = None):
    '''
    Parameters:
        peaks: all the peaks read in.
        num_poses: the number of positions on each strand that a read should occur.
        rpm: minimum number of reads.
        ratio: the ratio between the max(freads,rreads)/min(freads, rreads). Should be >= 1, thus only [1, ratio) will be accepted.
    Output:
        filtered_peaks: in the same format as the input peaks.
    '''
    filtered = {}
    for chrom in peaks:
        if chrom not in chroms:
            continue
        if chrom not in filtered:
            filtered[chrom] = []
        for p in peaks[chrom]:
            sum1 = p.data[0].sum()
            sum2 = p.data[1][:p.size].sum()
            if (p.data[0] > 0).sum() >= num_poses and (p.data[1][:p.size] > 0).sum() >= num_poses and sum1+sum2 >= rpm and not isReject(rejected, p.chrom, p.start + 1, p.end):
                if max(sum1, sum2) * 1.0 / min(sum1, sum2) < ratio:
                    filtered[chrom].append(p)

    return filtered

def stitchPeaks(peaks, stitchDist):
    '''
    Parameters:
        peaks : the peaks.
        stitchDist: the gap allowed between peaks.
    Output: a dictionary of the stitched peaks with chromosome as the key.
            {chrom1: [stitch1, stitched2, ...], chrom2:[]}
    '''
    stitched = {}
    index = 0 #track the id of the stitched peaks.
    for chrom in peaks:
        if len(peaks[chrom]) <= 0:
            continue
        if chrom not in stitched:
            stitched[chrom] = []
        lastEnd = peaks[chrom][0].end
        startIdx = 0
        for i in range(len(peaks[chrom])):
            if peaks[chrom][i].start >= lastEnd + stitchDist:
                stitched[chrom].append(StitchedPeaks(index, peaks[chrom][startIdx:i],stitchDist))
                index += 1
                startIdx = i
            lastEnd = peaks[chrom][i].end
        stitched[chrom].append(StitchedPeaks(index, peaks[chrom][startIdx:],stitchDist))
    return stitched



def loadBedGraph(reader, peaks):
    '''
    Parameters:
        reader: the BEDReader to read in bedgraph files, assuming the data is sorted.
        peaks: {chrom: [Peak, Peak...], chrom2: [Peak, Peak...]...}.
    Output: The data field in peaks will be filled.
        chroms: The list of chromosomes with data.
    '''
    chroms = []
    rn = reader.next
    def getNextLine():
        try:
            l = rn()
        except StopIteration:
            l = None
        return l
    line = getNextLine()
    while True:
        if line == None:
            break
        currChrom = line[0]
        if currChrom in peaks:
            chroms.append(currChrom)
            currPeaks = peaks[currChrom]
            i = 0
            tempLine = line
            while i < len(currPeaks) and tempLine != None and tempLine[0] == line[0]:
                p = currPeaks[i]
                p.prepareData()
                if p.end <= tempLine[1]:
                    i += 1
                elif tempLine[2] <= p.start:
                    tempLine = getNextLine()
                else:
                    p.addData(*tempLine[1:4])
                    if tempLine[2] < p.end:
                        tempLine = getNextLine()
                    elif p.end < tempLine[2]:
                        i += 1
                    else:
                        tempLine = getNextLine()
                        i += 1
            #initiate the remaining peaks
            while i < len(currPeaks):
                currPeaks[i].prepareData()
                i += 1
            while tempLine != None and tempLine[0] == line[0]:
                tempLine = getNextLine()
            line = tempLine
    return chroms


def pair( fp, rp, ulimit, dlimit):
    '''
    '''
    #Store the mates. (mate index, score)
    pairF = []
    pairR = []
    rstarts = []
    rends = []
    #preferences of the forward and reverse peaks.
    fprefer = []
    rprefer = []
    unpairedF = []
    defaultV = (-1,0)
    for f in fp:
        fprefer.append( [] )
        pairF.append(defaultV)
    for r in rp:
        rprefer.append( [] )
        pairR.append(defaultV)
        rstarts.append(r.start)
        rends.append(r.end)

    for i in range(len(fp)):
        currfp = fp[i]
        es = currfp.start - ulimit
        ee = currfp.end + dlimit

        si = bisect.bisect_left(rends, es) #find the left border of the reverse peaks to be included
        ei = bisect.bisect_right(rstarts, ee) #find the right border of the reverse peaks to survey

        for ri in range(si, ei):
            currrp = rp[ri]
            ccn = currfp.getCC( currrp )
            tempScore = max(ccn)
            fprefer[i].append( (ri, tempScore) )
            rprefer[ri].append( (i, tempScore) )
    for i in range(len(fp)):
        fprefer[i].sort(key=lambda k:(k[1],))
    for i in range(len(rp)):
        rprefer[i].sort(key=lambda k:(k[1],))
    unpairedF = range(len(fp))
    i = 0
    while len(unpairedF) > i:
        fi = unpairedF[i]
        #if fp[fi].start == 127571011:
        if len(fprefer[fi]) > 0:
            f = fprefer[fi][-1]
            ri = f[0]
            if pairR[ri][0] < 0 or f[1] > pairR[ri][1]:
                pairF[fi] = f
                # If the reverse peak is already paired, add the ex-forward peak to the unpairedF
                # And reset its pair information.
                if pairR[ri][0] > 0:
                    tempI = pairR[ri][0]
                    pairF[tempI] = defaultV
                    unpairedF.append(tempI)
                pairR[ri] = (fi, f[1])
            #fprefer[fi].remove(f)
            del fprefer[fi][-1]
        i += 1
    pairs = []
    singletons = []
    for i,p in enumerate(pairF):
        if p[0] >= 0:
            pairs.append(PeakPair(fp[i], rp[p[0]]))
        else:
            singletons.append(fp[i])
    for i,p in enumerate(pairR):
        if p[0] < 0:
            singletons.append(rp[i])

    return pairs, singletons

def writeData(stitched, peaks, prefix, total, ctotal):
    outStitched = open(prefix+'_stitched.xls','w')
    outStitched.write(getStitchHeader())
    outStitched.write('\t%s\t%s\n'%("Total_Reads", "Total_Ctrl",))
    for chrom in stitched:
        for s in stitched[chrom]:
            outStitched.write(str(s))
            outStitched.write('\t%d\t%d'%(total, ctotal,))
            outStitched.write('\n')
    outStitched.close()
    outPeaks = open(prefix+'_filtered.narrowPeak','w')
    for chrom in peaks:
        for p in peaks[chrom]:
            outPeaks.write(str(p))
            outPeaks.write('\n')
    outPeaks.close()

def getTotal(fwig):
    if fwig == None:
        return 0
    total = 0
    for chrom in fwig:
        total += fwig[chrom][:,1].sum()
    return total

def getOutput( options):
    output = options.output + '_n%d_m%.0f_r%.1f'%(options.nposes, options.rpm, options.ratio)
    return output

def plotDistHist(filtered, options):
    shifts = []
    oldShifts = []
    for chrom in filtered:
        for p in filtered[chrom]:
            p.getCC()
            p.getOldShift()
            shifts.append(p.shift)
            oldShifts.append(p.old_shift)

            if options.savefig >= 0:
                if options.rpm == options.savefig:
                    dirname = "%s/%d"%(getOutput(options),p.shift/10)
                    if not os.path.exists(dirname):
                        os.system("mkdir -p " + dirname)
                    p.plot(dirname)

    import pylab as pl
    pl.figure()
    pl.subplot(211)
    pl.hist(shifts, bins=100, range=(0,200))
    pl.subplot(212)
    pl.hist(oldShifts, bins=100, range=(0,200))
    pl.savefig(getOutput(options)+"_shifts.png",dpi=600)
    pl.cla()
    pl.clf()
    pl.close()


def run(options):
    #peaks = loadNarrowPeaks(BEDReader("RNAP-II_8WG16_SNU16_XO111_sh10_t40_peaks.narrowPeak",'narrowPeak'))
    peaks = loadNarrowPeaks(BEDReader(options.peakfile,'narrowPeak'))
    #fwig, rwig = addWigToPeak("/home/caofan/Downloads/MJF11_hg19/1_Bam/test_apex/RNAP-II_8WG16_SNU16_XO111_Forward.wig", "/home/caofan/Downloads/MJF11_hg19/1_Bam/test_apex/RNAP-II_8WG16_SNU16_XO111_Reverse.wig",peaks)
    fwig, rwig, chroms= addWigToPeak(options.fwig, options.rwig, peaks)
    #cfwig, crwig = addWigToPeak("/home/caofan/Downloads/MJF11_hg19/1_Bam/test_apex/SNU16_merged_Forward.wig", "/home/caofan/Downloads/MJF11_hg19/1_Bam/test_apex/SNU16_merged_Reverse.wig",peaks, True)
    cfwig = None
    crwig = None
    if options.cfwig and options.crwig:
        cfwig, crwig, chroms_c = addWigToPeak(options.cfwig, options.crwig, peaks, True)

    total = getTotal(fwig) + getTotal(rwig)
    ctotal = getTotal(cfwig) + getTotal(crwig)
    for r in [2,5,10,20]:
        options.rpm = r
        rpm = total*options.rpm/1000000
        print total, ' ', rpm, ' ',ctotal
        rejected = {}
        if options.rejectRegion is not None:
            rejected = loadRejectRegion(options.rejectRegion)
        filtered = filterPeaks(peaks, options.nposes, rpm, options.ratio, chroms, rejected)
        stitched = stitchPeaks(filtered, 12500)
        for chrom in stitched:
            for s in stitched[chrom]:
                s.getTotalCtrl(cfwig, crwig)
                s.getTotalSignal(fwig, rwig)
                s.getTotalSignalInConst()
                s.getTotalCtrlInConst()
        plotDistHist(filtered, options)
        writeData(stitched, filtered, getOutput(options), total, ctotal)
    import pylab as pl
    pl.clf()
    pl.close()




def getArgs():
    import argparse,sys
    parser = argparse.ArgumentParser(description="Stitch peaks.")
    parser.add_argument('peakfile',help="The peak file in narrowPeak format.")
    parser.add_argument('fwig', help="The forward strand wig file.")
    parser.add_argument('rwig', help="The reverse strand wig file.")
    parser.add_argument('--cfwig', help="The control forward strand wig file.")
    parser.add_argument('--crwig', help="The control reverse strand wig file.")
    parser.add_argument('-n','--nposes', type=int, default=0, help="The minimum number of locations a read could be found on each strand. Default: 0.")
    parser.add_argument('-m','--rpm', type=float, default=0, help="The minimum number of reads per million mapped reads a peak should have. Default: 0.")
    parser.add_argument('-r', '--ratio', type=float, default=sys.maxint,help="The maximum ratio between the forward and reverse strand reads. Default: %d"%(sys.maxint))
    parser.add_argument('-o', '--output', default="New", help="The prefix for output. If not specified, it will use 'New'.")
    parser.add_argument('-b', '--rejectRegion', help = "The blacklisted regions within which a peak should be rejected.")
    parser.add_argument('--savefig', type=int, default=-1, help="Whether to save the shift plots.")

    args = parser.parse_args()
    return args

if __name__=='__main__':
    options = getArgs()
    run(options)
