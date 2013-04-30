import os,sys
import pylab as plt
import argparse

COLORS = [(1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
        (0.0, 0.0, 0.0),
        (0.5882352941176471, 0.5098039215686274, 0.0),
        (1.0, 0.0, 1.0),
        (0.0, 1.0, 1.0),
        (0.39215686274509803, 0.39215686274509803, 0.39215686274509803)]


def loadData( filename, sep ):
    f = open( filename )
    data = []
    for r in f:
        if sep == None:
            data.append( [ float(i) for i in r.strip().split() ] )
        else:
            data.append( [ float(i) for i in r.strip().strip(sep).split(sep) ] )

    f.close()
    return data

def loadLabels( filename ):
    f = open( filename )
    labels = []
    for r in f:
        labels.append( r.strip() )
    f.close()
    return labels

def main( args ):
    data = loadData( args.input, args.sep )
    cdf, bins, patches = plt.hist( data, bins=args.bins, range=( args.min, args.max ), normed=1, cumulative=True )

    if args.label != None:
        labels = loadLabels( args.label )
    else:
        labels = [ str(i) for i in range( len( cdf ) )]
    

    centers = 0.5 * ( bins[1:] + bins[:-1] )
    plt.clf()
    for i,c in enumerate( cdf ):
        plt.plot( centers, c, color=COLORS[ i / args.perSample ], label = labels[i] + str(sum(data[i])))
    plt.legend()
    plt.title( args.title )
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser("Plot histogram for the data in a file.")
    parser.add_argument('-s', '--sep', help="The separater of the data.")
    parser.add_argument('-l', '--label', help="The file with the labels.")
    parser.add_argument('-b', '--bins', type=int, default=50, help="The number of bins")
    parser.add_argument('input', help="The file contains the data. Each row represent a sample.")
    parser.add_argument('min', type=int, help="Minimum for plot")
    parser.add_argument('max', type=int, help="Maximum for plot")
    parser.add_argument('-p', '--perSample', type=int, default=1, help="How many samples to share the same color.")
    parser.add_argument('-t', '--title', default= "", help="The title for the plot")

    args = parser.parse_args()
    main( args )

