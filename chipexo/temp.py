import test_pair_cc as tcc
import pylab as pl
def getOverlap(d1, d2, cc):
    overlap = []
    for i in range(len(cc)):
        overlap.append(min(i+1,d2.shape[0])-max(i+1-d1.shape[0],0))
    return overlap

def plotData(p, savefig=False):
    import numpy as np
    fdata = (p.fpeak.data - p.fpeak.data.mean())/p.fpeak.data.std()
    rdata = (p.rpeak.data - p.rpeak.data.mean())/p.rpeak.data.std()
    pl.subplot(231)
    pl.plot(range(p.fpeak.start, p.fpeak.end), fdata)
    pl.plot(range(p.rpeak.start, p.rpeak.end), rdata, color='r')
    pl.subplot(232)
    pl.plot(np.array(range(p.fpeak.start, p.fpeak.end)) + p.shift/2, fdata)
    pl.plot(np.array(range(p.rpeak.start, p.rpeak.end)) - p.shift + p.shift/2, rdata, color='r')

    overlap = getOverlap(p.fpeak.data, p.rpeak.data, p.ccn)
    pl.subplot(233)
    cc = np.zeros_like(p.ccn)
    for i in range(len(cc)):
        cc[i] = p.ccn[i] * overlap[i]
    shift2 = np.argmax(cc)+p.rpeak.start+1-p.fpeak.end
    pl.plot(np.array(range(p.fpeak.start, p.fpeak.end)) + shift2/2, fdata)
    pl.plot(np.array(range(p.rpeak.start, p.rpeak.end)) - shift2 + shift2/2, rdata, color='r')
    pl.subplot(234)
    pl.plot(overlap)
    pl.subplot(235)
    pl.plot(p.ccn)
    pl.subplot(236)
    pl.plot(cc)


    if savefig:
        pl.savefig('temppng_d_s/%s_%s_%d.png'%(p.fpeak.name, p.rpeak.name, p.shift))
        pl.clf()
    return p.shift, shift2

def loadAll():
    reader = tcc.BEDReader('../../macs_modified_out_control/MAX_sc-197_SNU16_XO111_+_peaks.narrowPeak', 'narrowPeak')
    fpeaks = tcc.loadNarrowPeaks(reader, 20)
    reader = tcc.BEDReader('../../macs_modified_out_control/MAX_sc-197_SNU16_XO111_-_peaks.narrowPeak', 'narrowPeak')
    rpeaks = tcc.loadNarrowPeaks(reader, 20)
    reader = tcc.BEDReader('../../macs_modified_out_control/MAX_sc-197_SNU16_XO111_+_ppois.bdg.bedGraph', 'bedGraph')
    tcc.loadBedGraph(reader, fpeaks)
    reader = tcc.BEDReader('../../macs_modified_out_control/MAX_sc-197_SNU16_XO111_-_ppois.bdg.bedGraph', 'bedGraph')
    tcc.loadBedGraph(reader, rpeaks)
    return fpeaks, rpeaks

def getShifts(pairs):
    shifts = []
    for p in pairs:
        shifts.append(p.shift)
    return shifts

def getCC(peak1, peak2):
    import numpy as np
    cc = np.correlate(data2, data1, 'full')
    for i in range(len(cc)):
        ccn[i] = cc[i] / (min(peak.start + i + 1, peak.end) - max(self.start + peak.start + i + 1 - self.end, peak.start))
    size1 = peak.end - peak.start
    size2 = self.end - self.start
    ratio = max(size1, size2) * 1.0 / min(size1, size2)
    ccn = ccn/ratio
    #x = np.argmax(ccn)
    #shift = peak.start + x + 1 - self.end
    return cc, ccn
