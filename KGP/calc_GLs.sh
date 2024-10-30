#!/bin/bash

module load bioinfo-tools 
module load samtools/1.5_debug

HREF=/proj/snic2020-2-10/1000AncientGenomes/ref_seqs/hs37d5.fa
CORES=8

ls /proj/snic2020-2-10/sllstore2017087/nobackup/private/reference_bias/refbias_GL/KGP/splitread_bams/*bam > tmp.bamfiles.txt

for pop in FIN YRI JPT
do
	for bamfile in $(grep $pop tmp.bamfiles.txt)
	do
		bambase=$(basename $bamfile .bam)
		python ../GL_calculator.py $bamfile snplist.txt counts.$bambase $HREF $pop.beagle $CORES
	done
	gzip -f $pop.beagle
	gzip -f $pop.beagle.corrected
done

