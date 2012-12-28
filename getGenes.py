import os,sys
sys.path.append('.')
import readBed
import findOverlap as fo

reader = readBed.BEDReader(sys.argv[1])
transcripts = []
for r in reader:
	transcripts.append(r)

transcripts.sort(key=lambda k:(k['name'],k['chrom'],int(k['chromStart']),int(k['chromEnd'])))

genes = {}
for t in transcripts:
	if t['name'] not in genes:
		genes[t['name']] = [[t['chrom'],t['chromStart'],t['chromEnd'],t['name'],t['score'],t['strand']],]
	else:
		if t['chrom'] != genes[t['name']][-1][0]:
			genes[t['name']].append([t['chrom'],t['chromStart'],t['chromEnd'],t['name'],t['score'],t['strand']])
		else:
			if fo.getOverlapOri(int(t['chromStart']),int(t['chromEnd']),int(genes[t['name']][-1][1]),int(genes[t['name']][-1][2])) > 0:
				if int(t['chromStart']) < int(genes[t['name']][-1][1]):
					genes[t['name']][-1][1] = t['chromStart']
				if int(t['chromEnd']) > int(genes[t['name']][-1][2]):
					genes[t['name']][-1][2] = t['chromEnd']
			else:
				genes[t['name']].append([t['chrom'],t['chromStart'],t['chromEnd'],t['name'],t['score'],t['strand']])

out = open(sys.argv[2],'w')
for g in genes:
	for r in genes[g]:
		out.write('\t'.join(r))
		out.write('\n')

out.close()
