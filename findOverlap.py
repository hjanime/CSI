import os,sys
import argparse
sys.path.append('.')
import readBed


def getOverlapOri(s1,e1,s2,e2):
	return max(0,min(e1,e2)-max(s1,s2))

def getOverlap(a,b,strand=False):
	if a['chrom'] != b['chrom'] or (strand and a['strand'] != b['strand']):
		return 0
	return max(0,min(int(a['chromEnd']),int(b['chromEnd'])) - max(int(a['chromStart']),int(b['chromStart'])))
	return getOverlap(int(a['chromStart']),int(a['chromEnd']),int(b['chromStart']),int(b['chromEnd']))

def getRecords(filename):
	f = readBed.BEDReader(filename)
	records = []
	for r in f:
		records.append(r)
	records.sort(key=lambda k:(k['chrom'],int(k['chromStart']),int(k['chromEnd'])))
	return records

def getArgs():
	parser = argparse.ArgumentParser(description="Find the overlapped features in two or more bed files.")
	parser.add_argument('infiles',metavar='I',type=str,nargs='+', help="Input prefixes. The bed files are prefix.bed and the coverage files are prefix.cov.txt")
	parser.add_argument('-f','--fraction',type=float,default=0,help="The overlapped region should have a minimum length of the specified fraction of the original feature")
	parser.add_argument('-r','--reciprocal',action='store_true',default=False,help="Whether to require the second feature also have the minimum length")
	parser.add_argument('-s','--strand',action='store_true',default=False,help="Whether to force strand when looking for overlappings")

	args = parser.parse_args()
	return args


def getNearFeatures(a,b,fraction,isReciprocal,useStrand):
	i = 0 #index for a
	j = 0 #index for b
	ao = {}
	bo = {}
	while i < len(a) and j < len(b):
		f1 = a[i]
		f2 = b[j]
		if f1['name'] not in ao:
			ao[f1['name']] = set()
		if f2['name'] not in bo:
			bo[f2['name']] = set()
		o = getOverlap(f1,f2,useStrand)
		#print f1['chromStart'],' ',f1['chromEnd'],' ',f2['chromStart'],' ',f2['chromEnd'],' ',o
		if o > 0 and o*1.0/(int(f1['chromEnd'])-int(f1['chromStart'])) > fraction:
			if (not isReciprocal) or o*1.0/(int(f2['chromEnd'])-int(f2['chromStart'])) > fraction:
				ao[f1['name']].add(f2['name'])
				bo[f2['name']].add(f1['name'])
		if (f1['chrom'] == f2['chrom'] and int(f1['chromEnd']) < int(f2['chromEnd'])) or f1['chrom'] < f2['chrom']:
			i += 1
		else:
			j += 1
	while i < len(a):
		if a[i]['name'] not in ao:
			ao[a[i]['name']] = set()
		i += 1
	while j < len(b):
		if b[j]['name'] not in bo:
			bo[b[j]['name']] = set()
		j += 1
	return ao,bo


def getNearMaxScore(a,b,covA,covB,fraction,isReciprocal,useStrand):
	i = 0 #index for a
	j = 0 #index for b
	ao = {}
	bo = {}
	while i < len(a) and j < len(b):
		f1 = a[i]
		f2 = b[j]
		if f1['name'] not in ao:
			ao[f1['name']] = [0,0]
		if f2['name'] not in bo:
			bo[f2['name']] = [0,0]
		o = getOverlap(f1,f2,useStrand)
		#print f1['chromStart'],' ',f1['chromEnd'],' ',f2['chromStart'],' ',f2['chromEnd'],' ',o
		if o > 0 and o*1.0/(int(f1['chromEnd'])-int(f1['chromStart'])) > fraction:
			if (not isReciprocal) or o*1.0/(int(f2['chromEnd'])-int(f2['chromStart'])) > fraction:
				if covB[f2['name']][0] > ao[f1['name']][0]:
					ao[f1['name']][0] = covB[f2['name']][0]
				if covB[f2['name']][1] > ao[f1['name']][1]:
					ao[f1['name']][1] = covB[f2['name']][1]
				#ao[f1['name']] = covB[f2['name']]
				if covA[f1['name']][0] > bo[f2['name']][0]:
					bo[f2['name']][0] = covA[f1['name']][0]
				if covA[f1['name']][1] > bo[f2['name']][1]:
					bo[f2['name']][1] = covA[f1['name']][1]
				#bo[f2['name']] = covA[f1['name']]
		if (f1['chrom'] == f2['chrom'] and int(f1['chromEnd']) < int(f2['chromEnd'])) or f1['chrom'] < f2['chrom']:
			i += 1
		else:
			j += 1
	while i < len(a):
		if a[i]['name'] not in ao:
			ao[a[i]['name']] = [0,0]
		i += 1
	while j < len(b):
		if b[j]['name'] not in bo:
			bo[b[j]['name']] = [0,0]
		j += 1
	return ao,bo


def getCoverage(filename):
	f = open(filename)
	coverage = {}
	for l in f:
		if l.startswith('#'):
			continue
		tokens = l.strip().split("\t")
		coverage[tokens[0]] = [int(i) for i in tokens[1:]]
	f.close()	
	return coverage

def getMaxCoverage(names,coverage):
	if len(names) == 0:
		return 0,0
	count = 0
	maxCov = 0
	for n in names:
		if coverage[n][0] > count:
			count = coverage[n][0]
		if coverage[n][1] > maxCov:
			maxCov = coverage[n][1]
	return count,maxCov
			
def generateOutput(records,overlaps,coverages,infiles,olType):
	for i in range(len(infiles)):
		generateSingleOutput(records,overlaps,coverages[i],i,infiles,olType)
			
def generateSingleOutput(records,overlaps,coverages,i,infiles,olType):
	prefix = infiles[i]
	ol = overlaps[i]
	out = open(prefix + ".overlap.out",'w')
	out.write("#chrom\tstart\tend\tname\tscore\tstrand\t%s_count"%(prefix,))
	for j in range(len(infiles)):
		if j != i:
			out.write("\t%s_count"%(infiles[j],))
	out.write("\t%s_maxCov"%(prefix,))
	for j in range(len(infiles)):
		if j != i:
			out.write("\t%s_maxCov"%(infiles[j],))
	out.write('\n')
	for r in records:
		out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t'%(r['chrom'],r['chromStart'],r['chromEnd'],r['name'],r['score'],r['strand'],coverages[r['name']][0]))
		maxCounts = []
		maxmaxCovs = []
		for j in range(len(infiles)):
			if j != i:
				tempOl = ol[j][r['name']]
				if olType == 'N':
					tempMaxCount,tempMaxmaxCov = getMaxCoverage(tempOl,coverages[j])
				else:
					tempMaxCount,tempMaxmaxCov = tempOl
				maxCounts.append(tempMaxCount)
				maxmaxCovs.append(tempMaxmaxCov)
		out.write('\t'.join([str(x) for x in maxCounts]))
		out.write("\t%s\t"%(coverages[r['name']][1],))
		out.write('\t'.join([str(x) for x in maxmaxCovs]))
		out.write('\n')

	out.close()

			

def main():
	args = getArgs()
	print args
	records = []
	overlaps = []
	coverages = []
	for filename in args.infiles:
	#	records.append(getRecords(filename + ".bed"))
		overlaps.append([])
	#	coverages.append(getCoverage(filename+'.cov.txt'))
	for i in range(len(args.infiles)):
		records1 = getRecords(args.infiles[i]+".bed")
		coverage1 = getCoverage(args.infiles[i] + ".cov.txt")
		for j in range(i,len(args.infiles)):
			print i,' ',j
			if i == j:
				overlaps[i].append({})
			else:
				records2 = getRecords(args.infiles[j]+".bed")
				coverage2 = getCoverage(args.infiles[j] + ".cov.txt")
				#ao,bo = getNearFeatures(records[i],records[j],args.fraction,args.reciprocal,args.strand)
				ao,bo = getNearMaxScore(records1,records2,coverage1,coverage2,args.fraction,args.reciprocal,args.strand)
				overlaps[i].append(ao)
				overlaps[j].append(bo)
		generateSingleOutput(records1,overlaps,coverage1,i,args.infiles,'S')
		overlaps[i] = []

	
if __name__ == '__main__':
	main()
