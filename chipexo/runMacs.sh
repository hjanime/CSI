large=20
small=1
tlarge=1
tsmall=40
root=~/Downloads/MJF11_hg19/1_Bam
for i in $*
do
    if [ ! -f $root/$i.unique.5.bed ]
    then 
        cat $root/$i.unique.+.5.bed  $root/$i.unique.-.5.bed > $root/$i.unique.5.bed
    fi
	macs2 callpeak -t ~/Downloads/MJF11_hg19/1_Bam/$i.unique.5.bed -f BED -g hs -n ${i} --bw 60 --verbose 2 --keep-dup all -q 0.001 --shiftsize ${large} --tsize ${tlarge}
	#macs2 callpeak -t ~/Downloads/MJF11_hg19/1_Bam/$i.unique.+.5.bed -f BED -g hs -n ${i}_+_sh${large}_t${tlarge} --bw 60 --verbose 2 --keep-dup all -q 0.001  --nomodel --shiftsize ${large} --tsize ${tlarge}
	#macs2 callpeak -t ~/Downloads/MJF11_hg19/1_Bam/$i.unique.+.5.bed -f BED -g hs -n ${i}_+_sh${small}_t${tsmall} --bw 60 --verbose 2 --keep-dup all -q 0.001  --nomodel --shiftsize ${small} --tsize ${tsmall}
	#macs2 callpeak -t ~/Downloads/MJF11_hg19/1_Bam/$i.unique.-.5.bed -f BED -g hs -n ${i}_-_sh${large}_t${tlarge} --bw 60 --verbose 2 --keep-dup all -q 0.001  --nomodel --shiftsize ${large} --tsize ${tlarge}
	#macs2 callpeak -t ~/Downloads/MJF11_hg19/1_Bam/$i.unique.-.5.bed -f BED -g hs -n ${i}_-_sh${small}_t${tsmall} --bw 60 --verbose 2 --keep-dup all -q 0.001  --nomodel --shiftsize ${small} --tsize ${tsmall}
	R CMD BATCH ${i}_model.r
done
