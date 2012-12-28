import os,sys
import argparse

def getMappedReads(filename):
	mapped = {}
	f = open(filename)
	for r in f:
		tokens = r.split()
		mapped[tokens[0].lower()] = long(tokens[1])
	f.close()
	return mapped

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Convert counts to RPKM values")
	parser.add_argument("infiles", metavar='I', type=str, nargs='+', help="Files that contain count values to be converted.")
	parser.add_argument("-m",'--mapped',type=str, help="File contains the counts of mapped reads.")
	parser.add_argument("-l",'--locFields',nargs='+',type=int, help="The fields that contains the location information. If 1 value is given, it will be treated as the length, otherwise please use Start followed by End field.")
	parser.add_argument('-f','--fields',type=int, nargs='+', help="The fields that has the values to be converted. Use 1-based coordinates.")
	args = parser.parse_args()

	mapCounts = getMappedReads(args.mapped)
	fields = args.fields
	print args.locFields
	for filename in args.infiles:
		print filename
		if 'rpkm' in filename:
			continue
		tag = filename.split('.')[0].lower()
		numReads = mapCounts[tag]
		f = open(filename)
		out = open(filename.replace(".bed",".rpkm.bed"),'w')
		for r in f:
			outStr = ''
			if r.startswith('#'):
				tokens  = r.split('\t')
				for i in fields:
					if i > 0 and i <= len(tokens):
						tokens[i-1] = 'rpkm'
				outStr = '\t'.join(tokens)

			elif r.startswith('track'):
				outStr = r
			else:
				tokens = r.strip().split('\t')
				length = 1
				if len(args.locFields) == 1:
					length = int(tokens[args.locFields[0]-1])
				else:
					length = int(tokens[args.locFields[1]-1]) - int(tokens[args.locFields[0]-1])

				for i in fields:
					if i>0 and i <= len(tokens):
						tokens[i-1] = str(int(tokens[i-1])*1.0*1e9/(length*numReads))
				outStr = '\t'.join(tokens) + "\n"
			out.write(outStr)
		out.close()
		
	
	
