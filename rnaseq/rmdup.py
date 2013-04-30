import os,sys

for bamfile in sys.argv[1:]:
	if bamfile.endswith("bam"):
		os.system("samtools rmdup %s %s"%(bamfile, bamfile.replace(".bam",".npcrd.bam")))


