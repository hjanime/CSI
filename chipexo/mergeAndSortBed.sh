for i in $*
do
    echo $i
    cat ${i}_+*peaks.out.bed ${i}_-*peaks.out.bed | sortBed > ${i}_peaks.bed
    cat ${i}_+*summits.out.bed ${i}_-*summits.out.bed | sortBed > ${i}_summits.bed
    
done
