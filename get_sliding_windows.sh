#!/usr/bin/env bash

# take human chr1 as an example
GENOME=./ref.fa
## make sliding windows: 500bp window size, 100bp step size
bedtools makewindows -b chidian.bed -w 500 -s 100 > Region.chr1.slide.chidian.bed
bedtools getfasta -tab -fi $GENOME -bed Region.chr1.slide.chidian.bed -fo Region.slide.chidian.input.txt
