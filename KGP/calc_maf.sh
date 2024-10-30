#!/bin/bash

module load bioinfo-tools
module load ANGSD/0.933


HREF=/proj/snic2020-2-10/1000AncientGenomes/ref_seqs/hs37d5.fa
ANC=/proj/snic2020-2-10/1000AncientGenomes/ref_seqs/hg19ancNoChr.fa

for pop in YRI FIN JPT
 do 
	angsd -doMaf 4 -beagle $pop.chip.selected.MAF0.2.nomiss.gl.gz -out $pop.true -fai $HREF.fai -anc $HREF
	 zcat $pop.beagle.gz | sed "s/chr//g" | gzip > beagle.tmp.gz
	zcat $pop.beagle.corrected.gz | sed "s/chr//g" | gzip > beagle.corrected.tmp.gz

	angsd -doMaf 4 -beagle beagle.tmp.gz -out $pop.def -fai $HREF.fai -anc $HREF
	angsd -doMaf 4 -beagle beagle.corrected.tmp.gz -out $pop.cor -fai $HREF.fai -anc $HREF
 done

