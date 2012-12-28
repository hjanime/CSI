import os,sys
sys.path.append('.')
import findOverlap as fo



def main():
	genes = fo.getRecords(sys.argv[2])
	annotatee = fo.getRecords(sys.argv[1])
	withinGeneType = sys.argv[3]
	withoutGeneType = sys.argv[4]
	fraction = float(sys.argv[5])
	ao,_ = fo.getNearFeatures(annotatee, genes,fraction,False,False)
	f = open(sys.argv[1])
	out = open(sys.argv[1]+'.annotation.txt','w')
	for r in f:
		if r.startswith('#'):
			out.write(r.strip() +"\ttype\tgene" + "\n")
		else:
			tokens = r.strip().split('\t')
			tempAnn = ao[tokens[3]]
			if len(tempAnn) > 0:
				annStr = withinGeneType + "\t" + ','.join(list(tempAnn))
			else:
				annStr = withoutGeneType
			tokens.append(annStr)
			out.write('\t'.join(tokens))
			out.write('\n')
	out.close()

	

if __name__=='__main__':
	main()
