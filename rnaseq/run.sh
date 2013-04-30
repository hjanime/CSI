GENES=~/Documents/genome/hg19/hg19.genes.bed
TRANSCRIPTS=~/Documents/genome/hg19/hg19.transcripts.bed
CDS=~/Documents/genome/hg19/hg19.cds.bed
ROOT=/media/HY-MJF/CaoFan/star_out
MAPPED=$ROOT/unique_mapped_counts.txt
for i in $*; do
	echo $i
	#samtools rmdup $i.bam $i.sorted.bam
	#samtools index $i.sorted.bam
	#python getUniqueHits.py $i.sorted.bam
	#sh getUniqueHits.sh $i
	python separateReads.py $TRANSCRIPTS $CDS $i.sorted.unique.bam
	#samtools index $i.sorted.unique.bam
	samtools index $i.sorted.unique.unknown.bam
	samtools index $i.sorted.unique.intragenic.bam
	samtools index $i.sorted.unique.exonic.bam
	samtools index $i.sorted.unique.intronic.bam
	python getClusters.py $i.sorted.unique.unknown.bam -o $i.unknown -l 76 -r 1 -N $MAPPED
	#python getClusters.py $i.sorted.unique.exonic.bam -o $i.exonic -l 74 -r 1
	#python getClusters.py $i.sorted.unique.intronic.bam -o $i.intronic -l 74 -r 1
	#python getCoverage.py $MAPPED $i
	#python merge.py $GENES $i 0.1
done

python findOverlap.py ../bam/Adult_rectum_fsorted.unknown ../bam/Adult_stomach_fsorted.unknown ../bam/Fetal_colon_fsorted.unknown ../bam/Fetal_heart_fsorted.unknown ../bam/Fetal_kidney_fsorted.unknown ../bam/Fetal_liver_fsorted.unknown ../bam/Fetal_lung_fsorted.unknown ../bam/Fetal_stomach_fsorted.unknown ../bam/SNU16_P13_fsorted.unknown ../bam/SNU16_P18_fsorted.unknown -f 0.1 
#python findOverlap.py Adult_stomach_fsorted.novel Fetal_stomach_fsorted.novel Fetal_colon_fsorted.novel -f 0.1

#for i in $*; do
#	python annotate.py $i.novel.overlap.out $GENES intron unknown 0.1
#	python annotate.py $i.exonic.overlap.out $GENES exon non_exon 0.1

#done


