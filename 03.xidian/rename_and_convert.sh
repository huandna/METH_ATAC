
for x in `ls ./*.sorted.bedGraph`; do
subanno_id=`echo ${x} | sed 's#.sorted.bedGraph##g' `
echo $subanno_id
bash ./replace_chr_in_bedgraph.sh ${subanno_id}.sorted.bedGraph  ${subanno_id}.rename.bedGraph chrs_paired.xls

LC_COLLATE=C sort -k1,1V -k2,2n -k3,3n sort -k1,1V -k2,2n -k3,3n ${subanno_id}.rename.bedGraph > ${subanno_id}.sorted.rename.bedGraph

bedGraphToBigWig ${subanno_id}.sorted.rename.bedGraph species_chrom_size.xls ${subanno_id}.bw
done



