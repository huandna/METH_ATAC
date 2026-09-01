
for x in `ls ./*.bedGraph`; do
subanno_id=`echo ${x} | sed 's#.bedGraph##g' `
echo $subanno_id
bedtools intersect -a $x -b chrom.size > tmp
sort -k1,1V -k2,2n -k3,3n tmp > ${subanno_id}.sorted.bedGraph
bedGraphToBigWig ${subanno_id}.sorted.bedGraph species_chrom_size.xls ${subanno_id}.bw
rm tmp;
done



