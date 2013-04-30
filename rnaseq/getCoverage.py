import os,sys
import pysam
import findOverlap as fo
import getClusters as gc
import count2rpkm as c2r
import re


def removeDuplicate(records):
	print 'start'
	records.sort(key=lambda k:(k['chrom'],int(k['chromStart']),int(k['chromEnd'])))
	currChrom = records[0]['chrom']
	currStart = records[0]['chromStart']
	currEnd = records[0]['chromEnd']
	currStrand = records[0]['strand']
	i = 1
	while i < len(records):
		if i%1000 == 0:
			print i
		if records[i]['chrom'] == currChrom and records[i]['chromStart'] == currStart and records[i]['chromEnd'] == currEnd and currStrand == records[i]['strand']:
			records[i]['chrom'] = 'd'
		else:
			currChrom = records[i]['chrom']
			currStart = records[i]['chromStart']
			currEnd = records[i]['chromEnd']
			currStrand = records[i]['strand']
		i += 1

def calculateCov(records, sam, outBed, outCov, numMapped):
	index = 1
	for r in records:
		if not re.match("^chr(\d+|X|Y|M)$",r['chrom']):
			continue
		count = 0
		for read in sam.fetch(r['chrom'],int(r['chromStart']),int(r['chromEnd'])-1):
			count += 1
		maxCov = 0
		for c in sam.pileup(r['chrom'],int(r['chromStart']), int(r['chromEnd']) -1):
			if c.n > maxCov:
				maxCov = c.n
		
		length = int(r['chromEnd']) - int(r['chromStart'])
		rpkm = gc.calRPKM(count,length,numMapped)
		outBed.write('\t'.join([r['chrom'],r['chromStart'],r['chromEnd'],str(index),str(rpkm),r['strand']]))
		outBed.write('\n')
		outCov.write('\t'.join([str(index),str(count),str(maxCov),str(rpkm)]))
		outCov.write('\n')
		index += 1

if __name__=="__main__":
	nameMap = {"exonic":"cds", "intronic":"noncds"}
	mapped = c2r.getMappedReads(sys.argv[1])
	for prefix in sys.argv[2:]:
		for k in nameMap:
			records = fo.getRecords("../genome/hg19/hg19."+nameMap[k]+".bed")
			removeDuplicate(records)
			sam = pysam.Samfile(prefix+".npcrd.unique."+k+".bam","rb")
			outBed = open(prefix+"."+k+".nonCluster.bed",'w')
			outCov = open(prefix + "."+k+".nonCluster.cov.txt",'w')
			calculateCov(records, sam, outBed, outCov, long(mapped[prefix.split('/')[-1].lower()]))
