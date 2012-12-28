import os,sys
sys.path.append('.')
import readBed
import findOverlap as fo

'''
Merges all the features.
'''
EXONIC_COLORS=["173,216,230","135,206,250","65,105,225","0,0,255","0,0,139"]
INTRONIC_COLORS=["220,220,220","192,192,192","128,128,128","80,80,80","0,0,0"]
UNKNOWN_COLORS=["255,255,154","238,220,130","205,155,29","139,139,0","139,107,0"]

def getColorId(cov):
	if cov > 100:
		return 4
	elif cov > 50:
		return 3
	elif cov > 30:
		return 2
	elif cov > 10:
		return 1
	else:
		return 0

def main():
	genes = fo.getRecords(sys.argv[1])
	fraction = float(sys.argv[-1])
	for prefix in sys.argv[2:-1]:
		out = open(prefix+'.all.bed','w')
		out.write("track name=\"%s\" visibility=2 itemRgb=\"On\"\n"%(prefix+'_all',))
		novel = fo.getRecords(prefix+".novel.bed")
		anns,_ = fo.getNearFeatures(novel,genes,fraction,False,False)
		f = open(prefix+'.novel.bed')
		for r in f:
			if r.startswith('#'):
				out.write(r.strip()+"\tthickStart\tthickEnd\titemRgb\n")
			else:
				tokens = r.strip().split('\t')
				tempAnn = anns[tokens[3]]
				tokens.append(tokens[1])
				tokens.append(tokens[2])
				cIndex = getColorId(int(tokens[4]))
				if len(tempAnn) > 0:
					tokens.append(INTRONIC_COLORS[cIndex])
				else:
					tokens.append(UNKNOWN_COLORS[cIndex])
				out.write('\t'.join(tokens))
				out.write('\n')

		f.close()
		exonic = fo.getRecords(prefix+".exonic.bed")
		anns,_ = fo.getNearFeatures(exonic,genes,fraction,False,False)
		f = open(prefix+".exonic.bed")
		for r in f:
			if r.startswith('#'):
				continue
			else:
				tokens = r.strip().split('\t')
				tempAnn = anns[tokens[3]]
				tokens.append(tokens[1])
				tokens.append(tokens[2])
				cIndex = getColorId(int(tokens[4]))
				tokens.append(EXONIC_COLORS[cIndex])
				tokens[3] = 'e'+tokens[3]
				out.write('\t'.join(tokens))
				out.write('\n')
		f.close()
		out.close()

if __name__ == '__main__':
	main()
				
					
					
		
		
