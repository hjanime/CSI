import os,sys
import pylab as pl
import argparse

if __name__=='__main__':
	parser = argparse.ArgumentParser("Generate box plots for the intrageneic features")
	parser.add_argument("infile",metavar='I',type=str,help="Input file of intragenic features annotations")
	parser.add_argument('title',metavar='T',type=str,help="Title for the box plot")
	args = parser.parse_args()

	f = open(args.infile)
	prefixes = f.readline().strip().split('\t')[5:]
	data = {}
	excounts = {}
	for p in prefixes:
		data[p] = []
		excounts[p] = 0
	for r in f:
		r = r.strip()
		if r == '':
			continue
		cov = r.split('\t')[5:]
		rname = r.split('\t')[4].lower()
		if 'rrna' in rname or 'trna' in rname:
			#print r
			continue
		if len(cov) < len(prefixes):
			#print r
			continue
		for i,p in enumerate(prefixes):
			if float(cov[i]) > 0.2:
				data[p].append(pl.log10(float(cov[i])))
				excounts[p] += 1
	fig = pl.figure()
	pl.boxplot([data[p] for p in prefixes])
	pl.ylim([-2,10])
	pl.ylabel('Log10 Expression Level (Log10(RPKM))')
	pl.xticks([i+1 for i in range(len(prefixes))], [i for i in prefixes],rotation=30)
	pl.title(args.title)
	fig.autofmt_xdate()
	pl.show()
