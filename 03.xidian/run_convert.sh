
for x in `ls ./*rename.bedGraph`; do
subanno_id=`echo ${x} | sed 's#.rename.bedGraph##g' `
echo $subanno_id
bedtools intersect -a $x -b chrom.size > tmp
LC_COLLATE=C sort -k1,1 -k2,2n tmp > ${subanno_id}.rename.sorted.bedGraph
bedGraphToBigWig ${subanno_id}.rename.sorted.bedGraph species_chrom_size.xls ${subanno_id}.bw
rm tmp;
done

