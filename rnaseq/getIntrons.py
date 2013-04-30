import os,sys
sys.path.append('.')
import findOverlap as fo

def getFeatures(filename):
	f = open(filename)
	features = []
	for r in f:
		if r.startswith("#"):
			continue
		tokens = r.strip().split('\t')
		tokens[1] = int(tokens[1])
		tokens[2] = int(tokens[2])
		features.append(tokens)
	features.sort(key=lambda k:(k[5],k[0],k[1],k[2])) # sort by strand first
	return features

def mergeFeatures(features):
	currChrom = features[0][0]
	currStart = features[0][1]
	currEnd = features[0][2]
	currStrand = features[0][5]
	merged= []
	count = 1
	for e in features:
		if e[5] == currStrand and currChrom == e[0] and fo.getOverlapOri(currStart,currEnd,e[1],e[2]) > 0:	
			if e[2] > currEnd:
				currEnd = e[2]
			if e[1] < currStart:
				currStart = e[1]
		else:
			merged.append([currChrom,currStart,currEnd,count,0,currStrand])
			currChrom = e[0]
			currStart = e[1]
			currEnd = e[2]
			currStrand = e[5]
			count += 1
	merged.append([currChrom,currStart,currEnd,count,0,currStrand])
	return merged

if __name__=="__main__":
	exons = mergeFeatures(getFeatures("../genome/hg19/hg19.cds.bed"))
	genes = mergeFeatures(getFeatures("../genome/hg19/hg19.transcripts.bed"))
	i = 0 #index for exons
	j = 0 #index for genes
	introns = []
	out = open("../genome/hg19/hg19.noncds.bed","w")
	count = 1
	tempStart = genes[0][1]
	while i < len(exons) and j < len(genes):
		if exons[i][5] < genes[j][5]:
			print 'strand i: ',i
			print exons[i]
			print genes[j]
			i += 1
		elif exons[i][5] > genes[j][5]:
			j += 1
			if j < len(genes):
				tempStart = genes[j][1]
		else:
			if exons[i][0] < genes[j][0] or (exons[i][0] == genes[j][0] and exons[i][1] < genes[j][1]):
				print 'chrom/loc i: ', i
				print exons[i]
				print genes[j]
				print genes[j-1]
				i += 1
			elif exons[i][0] > genes[j][0] or (exons[i][0] == genes[j][0] and exons[i][2] > genes[j][2]):
				if tempStart < genes[j][2]:
					tempIntron = [exons[i][0],tempStart,genes[j][2],count,0,genes[j][5]]
					#introns.append([exons[i][0],genes[j][1],exons[i][1],count,0,exons[i][5]])
					out.write('\t'.join([str(k) for k in tempIntron]) + "\n")
					count += 1
				j += 1
				if j < len(genes):
					tempStart = genes[j][1]
			else:
				if tempStart < exons[i][1]:
					tempIntron = [exons[i][0],tempStart,exons[i][1],count,0,exons[i][5]]
					#introns.append([exons[i][0],genes[j][1],exons[i][1],count,0,exons[i][5]])
					out.write('\t'.join([str(k) for k in tempIntron]) + "\n")
					count += 1
				tempStart = exons[i][2]
				if tempStart >= genes[j][2]:
					j += 1
					if j < len(genes):
						tempStart = genes[j][1]
				i += 1
	out.close()
