'''
Created on 12 Sep, 2013

@author: caofan
'''

import os,sys

sys.path.append('/home/caofan/Documents/scripts/chipexo/')

from exo.models import Peak

import wig
from dpgmm.dpgmm import DPGMM
from sklearn import mixture
import bisect
import numpy as np
import pylab as pl
import itertools

def getPeakWig(w, peak):
    if peak.chrom not in w:
        return None
    chromWig = w[peak.chrom]
    ws = bisect.bisect_left(chromWig[:,0],peak.start+1) #The start of the region
    we = bisect.bisect_right(chromWig[:,0], peak.end) #The right side of the regions
    axis_x = []
    x = [] #The data for processing
    for i in range(ws, we):
        axis_x.append(chromWig[i,0])
        x += [chromWig[i,0]] * int(chromWig[i,1])
    axis_x = np.array(axis_x)
    x = np.array(x)
    return axis_x, x, chromWig[ws:we,1]

if __name__ == '__main__':
    color_iter = itertools.cycle(['k','r','g','b','c','m','y'])
    FILTER_VAL = [99999,100000,1000000]
    fwig = wig.loadWig('/home/caofan/Downloads/MJF11_hg19/1_Bam/test_apex/MAX_sc-197_SNU16_XO111_Forward.wig', smooth=False)
    rwig = wig.loadWig('/home/caofan/Downloads/MJF11_hg19/1_Bam/test_apex/MAX_sc-197_SNU16_XO111_Reverse.wig', smooth=False)
    peaks = Peak.objects.filter(run=9).order_by('-size')
    
    for i in range(10):
        print peaks[i]
        for filter_val in FILTER_VAL:
            
            if peaks[i].strand == '+':
                axis_x, x, orig_y = getPeakWig(fwig, peaks[i])
            else:
                axis_x, x, orig_y = getPeakWig(rwig, peaks[i])
            model = DPGMM(1)
            skmodel = mixture.DPGMM(n_components=8,alpha=32,n_iter=10000)
            min_x = axis_x[0]
            axis_x = axis_x - min_x
            x = x-min_x
            data = []
            print min_x
            for v in x:
                model.add([v])
                data.append(v)
            #print model.data
            skmodel.fit(data)
            print "data mean: ", np.average(x)
            print "SK means: ",skmodel.means_, ' ', skmodel._get_covars(), ' ', skmodel.bic(x),' ', skmodel.converged_
            model.setPrior()
            print 'start'
            model.solveGrow()
            draw = model.intMixture()
            axis_y = []
            print 'solved'
            print 'v: ', model.v
            for idx,p in enumerate(draw[1]):
                temp = []
                for j in axis_x:
                    temp.append(p.prob([j]))
                #if p.getCovariance()[0]/draw[0][idx] <= filter_val:
                axis_y.append(temp)
                #print "mean, covariance:", p.mean, ' ', p.getCovariance(), ' ', draw[0][idx], ' ', p.getCovariance()[0]/draw[0][idx]
            pl.scatter(axis_x+min_x, orig_y)
            #print axis_x
            #print axis_x.shape
            #print (axis_x+min_x).shape
            for idx,(y,c) in enumerate(zip(axis_y,color_iter)):
                pl.scatter(axis_x+min_x, np.array(y)*len(x)*draw[0][idx], color=c)
            #pl.show()
            pl.savefig('%s_%d_%d_%s_%d.png'%(peaks[i].chrom, peaks[i].start, peaks[i].end, peaks[i].strand,filter_val))
            pl.clf()

            