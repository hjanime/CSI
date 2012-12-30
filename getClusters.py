import os,sys
import pysam
import argparse

RCUT = 0
RLENGTH = 0
USESTRAND = False
AINSERT = 0
MINSERT = 0
NUMMAPPED = 0

def getStartEnd(r):
	rStart = r.pos
	rEnd = r.pos + RLENGTH
	if r.mate_is_unmapped:
		return rStart,rEnd
	if r.is_reverse:
		if r.isize <= -MINSERT:
			rStart += -AINSERT + RLENGTH
		else:
			rStart += r.isize + RLENGTH
	if not r.is_reverse:
		if r.isize >0 and r.isize < MINSERT:
			rEnd = r.pos + r.isize
		else:
			rEnd = r.pos + AINSERT
	return rStart,rEnd

def printHeader(outBed,outCov):
	outBed.write("#Chrom\tstart\tend\tnumber\tscore\tstrand\n")
	outCov.write("#number\tcount\tmaxCov\tRPKM\n")


def printRecord(outBed,outCov,count,start,end,strand,chrom,index,maxCov):
	if count[strand] == -1 or count[strand] < RCUT:
		return index
	strandStr = '+'
	if strand == 1:
		strandStr = '-'
	index += 1
	rpkm = count[strand]*1.0*1e9/((end[strand]-start[strand])*NUMMAPPED)
	outBed.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(chrom,start[strand],end[strand],index,rpkm,strandStr))
	outCov.write("%s\t%s\t%s\t%s\n"%(index,count[strand],maxCov[strand],rpkm))
	count[strand] = -1
	start[strand] = -1
	end[strand] = -1
	return index

def getClusters(infile,outBed,outCov):
	if infile.endswith("bam"):
		f = pysam.Samfile(infile,"rb")
	else:
		f = pysam.Samfile(infile,"r")
	printHeader(outBed,outCov)
	curr_chrom = 0  #the chrom in pysam is represented by integers
	total = 0
	maxCov = [-1,-1]
	count = [-1,-1]
	coverage = [[0],[0]]
	start = [-1,-1]
	end = [-1,-1]
	actStart = [-1,-1]
	actEnd = [-1,-1] #Used for calculation of coverage.
 
	
	for r in f.fetch():
		curr_strand = 0
		if USESTRAND and r.is_reverse:
			curr_strand = 1
		rStart,rEnd = getStartEnd(r)

		if r.rname == curr_chrom:
			if (rStart < end[curr_strand] and rStart >= start[curr_strand]) or \
			   (rEnd < end[curr_strand] and rEnd > start[curr_strand]) or \
			   (rStart <= start[curr_strand] and rEnd >= end[curr_strand]):
				count[curr_strand] += 1
				
				if start[curr_strand] < 0:
					start[curr_strand] = rStart
					coverage[curr_strand] = [1 for i in range(RLENGTH)]
				else:
					diff = r.pos + RLENGTH - actEnd[curr_strand]
					tempCov = [1 for i in range(RLENGTH)]
					maxCov[curr_strand] = max(coverage[curr_strand])
					if diff <= RLENGTH:
						tempCov[:(RLENGTH-diff)] = [i+1 for i in coverage[curr_strand][diff:]]
					coverage[curr_strand] = tempCov
				
				if rEnd > end[curr_strand]:
					end[curr_strand] = rEnd
				if r.pos + RLENGTH > actEnd[curr_strand]:
					actEnd[curr_strand] = r.pos + RLENGTH
				if rStart < start[curr_strand]:
					start[curr_strand] = rStart
			else:
				maxCov[curr_strand] = max(maxCov[curr_strand],max(coverage[curr_strand]))
				total = printRecord(outBed,outCov,count,actStart,actEnd,curr_strand,f.getrname(curr_chrom),total,maxCov)
				count[curr_strand] = 1
				start[curr_strand] = rStart
				end[curr_strand] = rEnd
				actStart[curr_strand] = r.pos
				actEnd[curr_strand] = r.pos + RLENGTH
				maxCov[curr_strand] = 1
				coverage[curr_strand] = [1 for i in range(RLENGTH)]
		else:
			first = 0
			if start[1] > 0 and start[1] < start[0]:
				first = 1
			maxCov[0] = max(maxCov[0],max(coverage[0]))
			maxCov[1] = max(maxCov[1],max(coverage[1]))
			total = printRecord(outBed,outCov,count,actStart,actEnd,first,f.getrname(curr_chrom),total,maxCov)
			total = printRecord(outBed,outCov,count,actStart,actEnd,1-first,f.getrname(curr_chrom),total,maxCov)
			

			curr_chrom = r.rname
			curr_strand = 0
			if USESTRAND and r.is_reverse:
				curr_strand = 1

			count[curr_strand] = 1
			start[curr_strand] = rStart
			end[curr_strand] = rEnd
			actStart[curr_strand] = r.pos
			actEnd[curr_strand] = r.pos + RLENGTH
			coverage[curr_strand] = [1 for i in range(RLENGTH)]
			coverage[1-curr_strand] = [0 for i in range(RLENGTH)]

			
	first = 0
	if start[1] > 0 and start[1] < start[0]:
		first = 1
	total = printRecord(outBed,outCov,count,start,end,first,f.getrname(curr_chrom),total,maxCov)
	total = printRecord(outBed,outCov,count,start,end,1-first,f.getrname(curr_chrom),total,maxCov)
	f.close()
	outBed.close()
	outCov.close()


def main():
	parser = argparse.ArgumentParser(description="Get clusters from a sorted bam/sam file of mapped reads")
	parser.add_argument('infiles', metavar='I',type=str,nargs='+',help="Sorted sam/bam files to be processed.")
	parser.add_argument('-l','--readLength',type=int,required=True,help="The length of the reads")
	parser.add_argument('-o','--outfiles',type=str,nargs='*',help="Files to write output to. The sequence correspondes to the sequence of input files. If too few outputs, the rest will be output to stdout")
	parser.add_argument('-r','--rcut',type=int,default=0,help="Minimum support")
	parser.add_argument('-s','--strand',action='store_true',default=False,help='Whether to force strandnees')
	parser.add_argument('-i','--insert',type=int,default=150,help="The average insert size.")
	parser.add_argument('-m','--maxInsert',type=int,default=400,help='Maximum insert size.')
	parser.add_argument('-N','--numMapped',type=long,required=True,help='Number of unique reads mapped.')

	args = parser.parse_args()
	global RCUT 
	RCUT = args.rcut
	global RLENGTH
	RLENGTH = args.readLength
	global USESTRAND
	USESTRAND = args.strand
	global AINSERT
	AINSERT = args.insert
	global MINSERT
	MINSERT = args.maxInsert
	global NUMMAPPED
	NUMMAPPED = args.numMapped
	outCount = 0
	for infile in args.infiles:
		if outCount >= len(args.outfiles):
			out = sys.stdout
		else:
			outBed = open(args.outfiles[outCount]+'.bed','w')
			outCov = open(args.outfiles[outCount]+'.cov.txt','w')
			outCount += 1
		getClusters(infile,outBed, outCov)

if __name__=="__main__":
	main()
