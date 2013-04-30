GENES=../genome/hg19/hg19.genes.bed
TRANSCRIPTS=../genome/hg19/hg19.transcripts.bed
CDS=../genome/hg19/hg19.cds.bed
MAPPED=../unique_mapped_counts.txt
for i in $*; do
	echo $i
	#sh getUniqueHits.sh $i
	#samtools rmdup $i.bam $i.npcrd.bam
	#samtools index $i.npcrd.bam
	#python getUniqueHits.py $i.npcrd.bam
	#samtools index $i.npcrd.unique.bam
	#python separateReads.py $TRANSCRIPTS $CDS $i.npcrd.unique.bam
	#samtools index $i.npcrd.unique.unknown.bam
	#samtools index $i.npcrd.unique.intragenic.bam
	#samtools index $i.npcrd.unique.exonic.bam
	#samtools index $i.npcrd.unique.intronic.bam
	python getClusters.py $i.npcrd.unique.unknown.bam -o $i.unknown -l 76 -r 1 -N $MAPPED
	#python getClusters.py $i.npcrd.unique.exonic.bam -o $i.exonic -l 74 -r 1
	#python getClusters.py $i.npcrd.unique.intronic.bam -o $i.intronic -l 74 -r 1
	python getCoverage.py $MAPPED $i
	#python merge.py $GENES $i 0.1
done

#python findOverlap.py Adult_stomach_fsorted.exonic Fetal_stomach_fsorted.exonic Fetal_colon_fsorted.exonic -f 0.1 
#python findOverlap.py Adult_stomach_fsorted.novel Fetal_stomach_fsorted.novel Fetal_colon_fsorted.novel -f 0.1

#for i in $*; do
#	python annotate.py $i.novel.overlap.out $GENES intron unknown 0.1
#	python annotate.py $i.exonic.overlap.out $GENES exon non_exon 0.1

#done


