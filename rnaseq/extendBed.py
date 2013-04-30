import os,sys
import argparse
import readBed

if __name__=='__main__':
	parser = argparse.ArgumentParser("Extend the feature in a BED file")
	parser.add_argument('-e','--extend',type=int, default=10, help="The size to extend on each side")
	parser.add_argument('files',nargs='+',help="The input files")

	args= parser.parse_args()
	
	for infile in args.files:
		print infile
		f = readBed.BEDReader(infile)
		out = readBed.BEDWriter(infile.replace('.bed','_e%d.bed'%(args.extend,)))
		for r in f:
			r['chromStart'] = str(int(r['chromStart'])-args.extend)
			r['chromEnd'] = str(int(r['chromEnd'])+args.extend)
			out.writerow(r)
