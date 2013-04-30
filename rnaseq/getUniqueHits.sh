for i in $*; do
	echo $i
	echo "getUniqueHits"
	samtools view $i.npcrd.bam -bq 20 > $i.npcrd.unique.bam
	samtools index $i.npcrd.unique.bam
done
