
shiftsize=$1
tagsize=$2
mode=$3
OUTDIR=sh${shiftsize}_t${tagsize}_${mode}
SCRIPTDIR=~/Documents/scripts/chipexo

if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

python $SCRIPTDIR/getSummit.py -c --mode ${mode} --shift ${shiftsize} --tsize $tagsize --rpkm 10 --pileup 0.5 --half-ext -b ~/Documents/hg19.rejectRegion.bed -f fragsizes.txt -i ~/Downloads/MJF11_hg19/1_Bam/ -o $OUTDIR/ -w ~/Downloads/MJF11_hg19/1_Bam/test_apex/  MAX_sc-197_SNU16_XO111 MAX_sc-197_SNU16_XO211 RNAP-II_8WG16_HCT116_XO111 RNAP-II_8WG16_HCT116_XO211 RNAP-II_8WG16_SNU16_XO111 RNAP-II_8WG16_SNU16_XO21 -R ../final3/mappedCounts.txt

   python $SCRIPTDIR/mergeAndSort.py -i $OUTDIR/MAX_sc-197_SNU16_XO111_*_sh*_t*_peaks.out.xls -o $OUTDIR/MAX_sc-197_SNU16_XO111_peaks.out.xls -H
   python $SCRIPTDIR/mergeAndSort.py -i $OUTDIR/MAX_sc-197_SNU16_XO211_*_sh*_t*_peaks.out.xls -o $OUTDIR/MAX_sc-197_SNU16_XO211_peaks.out.xls -H
   python $SCRIPTDIR/mergeAndSort.py -i $OUTDIR/RNAP-II_8WG16_HCT116_XO111_*_sh*_t*_peaks.out.xls -o $OUTDIR/RNAP-II_8WG16_HCT116_XO111_peaks.out.xls -H
   python $SCRIPTDIR/mergeAndSort.py -i $OUTDIR/RNAP-II_8WG16_HCT116_XO211_*_sh*_t*_peaks.out.xls -o $OUTDIR/RNAP-II_8WG16_HCT116_XO211_peaks.out.xls -H
   python $SCRIPTDIR/mergeAndSort.py -i $OUTDIR/RNAP-II_8WG16_SNU16_XO111_*_sh*_t*_peaks.out.xls -o $OUTDIR/RNAP-II_8WG16_SNU16_XO111_peaks.out.xls -H
   python $SCRIPTDIR/mergeAndSort.py -i $OUTDIR/RNAP-II_8WG16_SNU16_XO21_*_sh*_t*_peaks.out.xls -o $OUTDIR/RNAP-II_8WG16_SNU16_XO21_peaks.out.xls -H

  sh mergeAndSortBed.sh $OUTDIR/MAX_sc-197_SNU16_XO111 $OUTDIR/MAX_sc-197_SNU16_XO211 $OUTDIR/RNAP-II_8WG16_HCT116_XO111 $OUTDIR/RNAP-II_8WG16_HCT116_XO211 $OUTDIR/RNAP-II_8WG16_SNU16_XO111 $OUTDIR/RNAP-II_8WG16_SNU16_XO21 

  java -jar ~/Documents/annotate.jar -a ~/hg19_wgGencodeAttrsV14.txt -f MACS -r ~/hg19_wgEncodeGencodeCompV14.txt -i ~/Desktop/macs_out/final6/${OUTDIR}/*1_peaks.out.xls -b ~/Documents/hg19.rejectRegion.bed
  java -jar ~/Documents/annotate.jar -a ~/hg19_wgGencodeAttrsV14.txt -f MACSSUMMIT -r ~/hg19_wgEncodeGencodeCompV14.txt -i ~/Desktop/macs_out/final6/${OUTDIR}/*1_peaks.out.xls -b ~/Documents/hg19.rejectRegion.bed

   for j in MAX_sc-197_SNU16_XO111 MAX_sc-197_SNU16_XO211 RNAP-II_8WG16_HCT116_XO111 RNAP-II_8WG16_HCT116_XO211 RNAP-II_8WG16_SNU16_XO111 RNAP-II_8WG16_SNU16_XO21
   do
       echo "pairing $OUTDIR/${i}"
       python $SCRIPTDIR/pair_shape.py $OUTDIR/${j}_peaks.bed ~/Downloads/MJF11_hg19/1_Bam/test_apex/${j}_Forward.wig  ~/Downloads/MJF11_hg19/1_Bam/test_apex/${j}_Reverse.wig $OUTDIR/${j}
   done


#echo "cp _peaks.bed /media/HY-MJF/chip-exo/final/annotation/peaks/"
#cp *1_peaks.bed /media/HY-MJF/chip-exo/final/bedfiles/peaks/
#echo "cp _summits.bed /media/HY-MJF/chip-exo/final/anntation/summits/"
#cp *1_summits.bed /media/HY-MJF/chip-exo/final/bedfiles/summits/

#cp *1_peaks.out.xls.annotation.macssummit /media/HY-MJF/chip-exo/final/annotation/summits/
#cp *1_peaks.out.xls.annotation.macs /media/HY-MJF/chip-exo/final/annotation/peaks/
#cp *XO111_peaks.out.xls *XO211_peaks.out.xls *XO21_peaks.out.xls *.out.bed /media/HY-MJF/chip-exo/final/

