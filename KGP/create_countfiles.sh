#!/bin/bash


module load bioinfo-tools
module load samtools/1.5_debug
module load bwa/0.7.17
module load python/3.7.2
module load pysam/0.15.3-python3.7.2
module load biopython/1.76-py3


BAMFILE=$1
CORES=4
REFGENOME=/proj/snic2020-2-10/1000AncientGenomes/ref_seqs/hs37d5.fa
TPED=pseudohaploid.all.tped

bambase=$(basename $BAMFILE .bam)

samtools view -L snplist.txt.bed -h -b -u $BAMFILE > $bambase.bedfilter.bam

samtools index	$bambase.bedfilter.bam

samtools calmd $bambase.bedfilter.bam $REFGENOME | ../modify_read_alternative.py snplist.txt | samtools view -b -u - > $bambase.modified.bam

bwa aln -l 16500 -n 0.01 -o 2 -t $CORES $REFGENOME -b $bambase.modified.bam \
| bwa samse $REFGENOME - $bambase.modified.bam  \
| samtools view -F 4 -h - \
| samtools sort -o $bambase.modified.remap.bam -

samtools index $bambase.modified.remap.bam

samtools calmd -b $bambase.modified.remap.bam $REFGENOME > $bambase.modified.remap.MD.bam
samtools index $bambase.modified.remap.MD.bam

samtools merge $bambase.merged.bam $bambase.modified.remap.MD.bam $bambase.bedfilter.bam
samtools index $bambase.merged.bam




python ../count_reads_mpB.py $bambase.merged.bam $REFGENOME $TPED snplist.txt.bed > counts.$bambase

rm $bambase.merged.bam $bambase.modified.remap.MD.bam $bambase.modified.remap.bam $bambase.modified.bam $bambase.bedfilter.bam
rm ${bambase}*bai
