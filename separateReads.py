import os,sys

if __name__=="__main__":
	for filename in sys.argv[2:]:
		os.system("intersectBed -abam %s -b %s -v > %s"%(filename.replace(".bam",".unique.bam"), sys.argv[1], filename.replace(".bam",".unique.unknown.bam")))
		os.system("intersectBed -abam %s -b %s -u > %s"%(filename.replace(".bam",".unique.bam"), sys.argv[1], filename.replace(".bam",".unique.intragenic.bam")))
		os.system("intersectBed -abam %s -b %s -u > %s"%(filename.replace(".bam",".unique.intragenic.bam"), sys.argv[1], filename.replace(".bam",".unique.exonic.bam")))
		os.system("intersectBed -abam %s -b %s -v > %s"%(filename.replace(".bam",".unique.intragenic.bam"), sys.argv[1], filename.replace(".bam",".unique.intronic.bam")))

	
