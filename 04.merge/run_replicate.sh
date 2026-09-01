#!/usr/bin/sh
python ./merge_methyl_replicates.py \
	-i ./methyl_bw -c /home/yaoxinw/workdir/project/11.nvwa/01.chidian/species_chrom_size.xls \
	-o ./merged_methyl_uniq \
	-t 0.05 \
	-b 10 \
	-n 24 > log_rep 2>&1

