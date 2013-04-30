import os,sys
import argparse



if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Convert GEM event predictions to gff format")
	parser.add_argument('infiles', metavar='I',nargs='+',help="Input files.")
	parser.add_argument('-e','--extend',type=int,default=0,help="To extend the feature to the left and right by the specified amount")
	args = parser.parse_args()
	def gem2gff(filename):

		count = 1
		f = open(filename)
		out = open(filename[:filename.rindex('.')]+'.gff','w')
		header = f.readline()
		names = header.strip().split('\t')
		for r in f:
			outArray = []
			tokens = r.strip().split('\t')
			tokens = [t.strip() for t in tokens]
			chrom,pos = tokens[0].split(':')
			chrom = 'chr'+chrom.upper()
			pos = int(pos)
			end = pos + 1
			pos -= args.extend
			end += args.extend
			score = tokens[1]
			description = 'Control=%s; Fold=%s; Expected=%s; Q_-lg10=%s; P_-lg10=%s; P_poiss=%s; IPvsEMP=%s; IPvsCTR=%s; KmerGroup=%s; KG_hgp=%s; Strand=%s'%(tokens[2],tokens[3],tokens[4],tokens[5],tokens[6],tokens[7],tokens[8],tokens[9],tokens[10],tokens[11],tokens[12])
			outArray += [chrom,'gem','.',str(pos),str(end),score,'.','.',description]
			out.write('\t'.join(outArray))
			out.write('\n')
		f.close()
		out.close()

		
	for f in args.infiles:
		gem2gff(f)
