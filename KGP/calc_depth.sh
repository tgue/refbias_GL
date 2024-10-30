#!/bin/bash

module load bioinfo-tools
module load bcftools/1.5 
module load samtools/1.5_debug
module load BEDTools


REFGENOME=/proj/snic2020-2-10/1000AncientGenomes/ref_seqs/hs37d5.fa
bamfile=$1

#cat bamlist.txt | while read bamfile; do 
samtools view -b -q 20 $bamfile | genomeCoverageBed -ibam - -g $REFGENOME > $TMPDIR/tmp.cov
coverage=$(grep genome $TMPDIR/tmp.cov | awk '{NUM+=$2*$3; DEN+=$3} END {print NUM/DEN}')
coverage_mt=$(grep MT $TMPDIR/tmp.cov | awk '{NUM+=$2*$3; DEN+=$3} END {print NUM/DEN}')
samtools view -q 20 $bamfile | awk '{print length($10)}' > $TMPDIR/length.txt
length=$(awk -v column_number=1 '{ sum += $column_number; count++ } END { printf sum / count }' $TMPDIR/length.txt)
min=$(sort -n $TMPDIR/length.txt | head -n 1)
max=$(sort -n $TMPDIR/length.txt | tail -n 1)

echo $bamfile $coverage $coverage_mt $length $min $max >> depth_stats.txt
#done

