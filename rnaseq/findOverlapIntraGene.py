import os,sys
import readBed

def printHeader(out,args):
	out.write("chrom\tchromStart\tchromEnd\tstrand\tname")
	for a in args:
		out.write("\t"+a.split('/')[-1])
	out.write('\n')

def getKey(r):
	return '\t'.join([r['chrom'],r['chromStart'],r['chromEnd'],r['strand']])

def init(reader):
	result = {}
	for r in reader:
		if r['chrom'] == 'chr1' and r['chromStart'] == 11873:
			print r
		tempkey = getKey(r)
		result[tempkey] = [':'.join(r['name'].split('///')[0:2]),]
	return result

def main():
	ereader = readBed.BEDReader('../genome/hg19/hg19.cds.bed')
	ireader = readBed.BEDReader('../genome/hg19/hg19.noncds.bed')
	exons = init(ereader)
	introns = init(ireader)
	for prefix in sys.argv[1:]:
		print prefix
	for prefix in sys.argv[1:]:
		ereader = readBed.BEDReader(prefix+".exonic.nonCluster.bed")
		for r in ereader:
			tempkey = getKey(r)
			exons[tempkey].append(r['score'])
		ireader = readBed.BEDReader(prefix+".intronic.nonCluster.bed")
		for r in ireader:
			tempkey = getKey(r)
			introns[tempkey].append(r['score'])
	out = open("exons.overlap.out",'w')
	printHeader(out,sys.argv[1:])
	for e in exons:
		out.write('\t'.join([e,]+exons[e]))
		out.write('\n')
	out.close()
	out = open('introns.overlap.out','w')
	printHeader(out,sys.argv[1:])
	for i in introns:
		out.write('\t'.join([i,]+introns[i]))
		out.write('\n')
	out.close()


if __name__=='__main__':
	main()
