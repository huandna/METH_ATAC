#!/usr/bin/env bash
#conda activate prokka_env

GENOME=./ref.fa
python ../get_genome.bed.py --input $GENOME --output ./andai.bed
## make sliding windows: 500bp window size, 100bp step size
bedtools makewindows -b andai.bed -w 500 -s 100 > Region.chr1.slide.andai.bed
bedtools getfasta -tab -fi $GENOME -bed Region.chr1.slide.andai.bed -fo Region.slide.andai.input.txt
python ../pred_Human_chr.py --region ./Region.slide.andai.input.txt --anno ../Female_mus_anno_5w.txt


# get bedgraph
for x in `ls ./Predict_track/merged_*`; do 
cat $x | perl -ne '@a=($_ =~ /(.*):(\d+)-(\d+)\t(.*)/); print "".join("\t", $a[0],  $a[1]+250-50, $a[1]+250+50, $a[3]), "\n";' > `echo ${x}|sed 's/\.txt//g'`.bedGraph;
done

#get bw
ln -s ./andai.bed chrom.size
awk '{print $1"\t"$3}' andai.bed > species_chrom_size.xls
species_chrom_size=species_chrom_size.xls
#cat ${species_chrom_size} |awk '{print $1"\t"0"\t"$2}' > chrom.size


for x in `ls ./Predict_track/merged_*.bedGraph`; do 
subanno_bedgraph=`echo ${x} | sed 's#merged_##g' `
subanno_id=`echo ${subanno_bedgraph} | sed 's#.bedGraph##g' `
echo $subanno_id
bedtools intersect -a $x -b chrom.size >tmp
sort -k1,1V -k2,2n -k3,3n tmp >$subanno_bedgraph 
bedGraphToBigWig $subanno_bedgraph ${species_chrom_size} ${subanno_id}.bw
rm tmp;
done

mkdir bw
mv ./Predict_track/*.bw ./bw/
#rm -rf Predict_track/
