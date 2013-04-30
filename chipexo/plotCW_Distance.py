import os,sys
import numpy as np
import pylab as plt

COLORS = [(1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
        (0.0, 0.0, 0.0),
        (0.5882352941176471, 0.5098039215686274, 0.0),
        (1.0, 0.0, 1.0),
        (0.0, 1.0, 1.0),
        (0.39215686274509803, 0.39215686274509803, 0.39215686274509803)]


def plotForAll( dists, labels ):
    print len(dists)
    pdf, bins, patches = plt.hist( dists, bins = 50, range=(0,100), normed=1)
    centers = 0.5 * ( bins[1:] + bins[:-1] )
    print len(centers)
    plt.clf()
    for i,c in enumerate( pdf ):
        plt.plot( centers, c, color=COLORS[ i % len(COLORS) ], label = labels[i] )
    plt.legend()
    plt.title("CW distance")
    plt.show()

def main():
    dists = []
    fns = []
    for fn in sys.argv[1:]:
        f = open( fn )
        tempDists = []
        for r in f:
            if r.startswith('#'):
                continue
            else:
                tempDists.append( float( r.strip().split('=')[-1] ) )
        f.close()
        dists.append(tempDists)
        fns.append( fn )
    plotForAll( dists, fns )


if __name__=="__main__":
    main()
