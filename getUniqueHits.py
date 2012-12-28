import os,sys
import pysam


def main():
	for filename in sys.argv[1:]:
		samfile = pysam.Samfile(filename,"rb")
		uniquefile = pysam.Samfile(filename.replace(".bam",".unique.bam"),"wb",template=samfile)
		total = 0
		uniqueCount = 0
		for read in samfile.fetch():
			total += 1
			if not read.is_unmapped:
				for t in read.tags:
					if t[0] == 'XT' and (t[1] == 'U' or t[1] == 'M'):
						uniqueCount += 1
						uniquefile.write(read)
		samfile.close()
		uniquefile.close()
		print "Unique/Total"
		print uniqueCount,'/',total
		os.system("intersectBed -abam %s -b %s -v > %s"%(filename.replace(".bam",".unique.bam"), "genome/hg19/hg19.exons.bed", filename.replace(".bam",".unique.novel.bam")))
		os.system("intersectBed -abam %s -b %s -u > %s"%(filename.replace(".bam",".unique.bam"), "genome/hg19/hg19.exons.bed", filename.replace(".bam",".unique.exonic.bam")))

main()
