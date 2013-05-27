import os,sys
import numpy as np
import re
import pylab as pl
from matplotlib.ticker import MultipleLocator
import argparse
import bisect

CATEGORIES = ['splicing','intronic','utr3','utr5','exonic','intergenic','upstream','downstream']
CATEGORIES_S = ['genic','distal']
COLORS = [(255,0,0),(230,127,0),(0,0,200),(200,200,0),(105,0,180),(0,180,0),(193,0,255),(140,140,140)]
COLORS_S = [(255,127,0),(140,140,140)]
MAP = {'splicing': 0, 'intronic': 1, 'utr3':2, 'utr5':3,'exonic':4, 'intergenic':5, 'upstream':6, 'downstream':7}
MAP_S = {'upstream': 0,'downstream': 0,'intronic': 0, 'utr3':0, 'utr5':0,'splicing': 0, 'exonic':0, 'intergenic':1, }
NEARTHR = 10000

def normColor(colors):
    result = []
    for c in colors:
        r = c[0]*1.0/255
        g = c[1]*1.0/255
        b = c[2]*1.0/255
        result.append((r,g,b))
    return result

def getTitle(filename):
    tokens = filename.split('/')[-1].split('_')
    if len(tokens) < 12:
        return filename.split('/')[-1].split('.')[0]
    used = tokens[1:4]
    used.append(tokens[7])
    used.append(tokens[10])
    used.append(tokens[11])
    return '_'.join(used)


def plotFigure(filename,data,cates,colors,xlabel,typefig,nbins=200):
    #out = pl.hist(data,nbins,normed=1,range=(0,100),histtype='bar',label=cates,color=colors,edgecolor='none')
    y = []
    binEdges = []
    total = []
    for d in data:
        total.append(len(d))
        tempy,binEdges = np.histogram(d,bins=range(201),normed=1)
        y.append(tempy)
    
    #print y
    #print binEdges
    binCenters = 0.5*(binEdges[1:]+binEdges[:-1])
    fig = pl.figure()
    for i,tempy  in enumerate(y):
        pl.plot(binCenters,tempy,'-', color=colors[i],label=cates[i]+" (%d)"%(total[i],))
    ax=fig.add_subplot(111)
    pl.grid(b=True,which='major', color='g', linestyle='--')
    pl.ylabel('Frequency',fontsize=16)
    pl.xlabel(xlabel,fontsize=16)
    ax.xaxis.set_major_locator(MultipleLocator(20.0))
    ax.xaxis.set_minor_locator(MultipleLocator(1.0))
    pl.legend()
    pl.title(getTitle(filename) + " " + xlabel, fontsize=10)
    pl.show()
    #pl.savefig(filename.split('.')[0]+"."+typefig+".png",dpi=600)
    

def plotBox(filename, data, intervals):
    xlabels = []
    for i in range( len(intervals) - 1 ):
        if i != len(intervals ) - 2:
            xlabels.append( '[%d,%d) %d'%(intervals[i], intervals[i+1],len(data[i])))
        else:
            xlabels.append( '[%d,%d] %d'%(intervals[i], intervals[i+1],len(data[i])))

    pl.boxplot(data) 
    pl.xticks( [i+1 for i in range( len( intervals) - 1 )], xlabels)
    pl.show()


def main( args ):
    colors = normColor(COLORS)
    colors_s = normColor(COLORS_S)
    for filename in args.inputs:
        print filename
        title = getTitle(filename)
        #exp = []
        dist = []
        #exp_s = []
        dist_s = []

        for c in CATEGORIES:
            #exp.append([])
            dist.append([])
        for c in CATEGORIES_S:
            #exp_s.append([])
            dist_s.append([])
        f = open(filename)
        if args.header:
            f.readline()
        lines = []
        for r in f:
            lines.append(r)
        f.close()
        #lines.sort(key=lambda k:(-float(k.split('\t')[7])))
        scores = []
        #remove duplicates from the interval list.
        args.intervals = list( set( args.intervals ) )
        if args.scoreCol and len( args.intervals ) > 1:
            args.intervals.sort()
            for i in range( len( args.intervals ) - 1 ):
                scores.append([])
            print args.intervals
        for r in lines:
            tokens = r.strip().split('\t')
            #read_count = float(tokens[])
            pair_dist = int(tokens[ args.distCol ].split('=')[-1])
            cates = set(tokens[ args.typeCol ].lower().split(':'))
            for c in cates:
                c = c.replace('ncrna_','') #remove the ncrna tag

                dist[ MAP[c] ].append( pair_dist )
                dist_s[ MAP_S[c] ].append( pair_dist )
            if args.scoreCol and len( args.intervals ) > 1:
                pair_score = float( tokens[ args.scoreCol ] )
                insert_pos_l = bisect.bisect_left( args.intervals, pair_dist )
                insert_pos_r = bisect.bisect_right( args.intervals, pair_dist )
                if pair_dist >= args.intervals[0] and pair_dist <= args.intervals[-1]:
                    if pair_dist == args.intervals[0]:
                        scores[0].append( pair_score )
                    else:
                        scores[ insert_pos_l -1 ].append( pair_score )


        f.close()
        
        #plotFigure(filename, exp, CATEGORIES, colors, 'Read Count', 'exp')
        
        
        plotFigure(filename, dist[1:], CATEGORIES[1:], colors, 'Distance (bp)', 'dist', 50)
        

        plotFigure(filename, dist_s, CATEGORIES_S, colors_s, 'Distance (bp)', 'dist_s', 50)

        if len( scores ) > 0:
            plotBox(filename, scores, args.intervals)


        #plotFigure(filename, exp_s, CATEGORIES_S, colors_s, 'Read Count', 'exp_s')



if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Plot the peak pair distances by categories.")
    parser.add_argument("inputs", nargs='+', help="The input files.")
    parser.add_argument('--typeCol', required=True, type=int, help="The index of the column for the annotation type. Index start from 0.")
    parser.add_argument('--scoreCol', type=int, help="The index of the column for the score. Index start from 0.")
    parser.add_argument('--intervals', type=int, nargs='+', help="The edges of the intervals for the box plots of scores, separated by space. For example: 0 10 20 30 will give interval [0,10), [10,20),[20,30].")
    parser.add_argument('--distCol', required=True, type=int, help="The index of the column for the distance. Index start from 0.")
    parser.add_argument('--header', action='store_true', default=False, help="Whether the first line is header.")
    args = parser.parse_args()
    main( args)

