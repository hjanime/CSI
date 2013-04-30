import os,sys

if __name__=='__main__':
	for filename in sys.argv[1:]:
		print filename
		f = open(filename)
		outWriter = open(filename.replace('.gff','.human'),'w')
		for r in f:
			tokens = r.strip().split('\t')
			out = []
			out.append(tokens[0][3:])
			out.append(tokens[3])
			out.append(tokens[4])
			out.append('0')
			out.append('0')
			out.append(tokens[5])
			out.append(tokens[-1])
			outWriter.write('\t'.join(out)+'\n')
		outWriter.close()
		
