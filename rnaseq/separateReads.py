import os,sys

if __name__=="__main__":
	for filename in sys.argv[3:]:
		os.system("intersectBed -abam %s -b %s -v > %s"%(filename, sys.argv[1], filename.replace(".bam",".unknown.bam")))
		os.system("intersectBed -abam %s -b %s -u > %s"%(filename, sys.argv[1], filename.replace(".bam",".intragenic.bam")))
		os.system("intersectBed -abam %s -b %s -u > %s"%(filename.replace(".bam",".intragenic.bam"), sys.argv[2], filename.replace(".bam",".exonic.bam")))
		os.system("intersectBed -abam %s -b %s -v > %s"%(filename.replace(".bam",".intragenic.bam"), sys.argv[2], filename.replace(".bam",".intronic.bam")))

	
