#!/bin/bash

module load bioinfo-tools
module load ANGSD/0.933



HREF=/proj/snic2020-2-10/1000AncientGenomes/ref_seqs/hs37d5.fa
ANC=/proj/snic2020-2-10/1000AncientGenomes/ref_seqs/hg19ancNoChr.fa

for subgenome in European African EastAsian
do
	for pop in FIN YRI JPT
	do
		angsd -doMaf 4 -beagle $pop.chip.selected.MAF0.2.nomiss.gl.gz -out $pop.true.$subgenome -fai $HREF.fai -anc $HREF -sites anc.hg19.$subgenome.bed
		angsd -doMaf 4 -beagle $pop.pseudohaploid.gl.gz -out $pop.ph.$subgenome -fai $HREF.fai -anc $HREF -sites anc.hg19.$subgenome.bed
		zcat $pop.beagle.gz | sed "s/chr//g" | gzip > $subgenome.tmp.gz
		angsd -doMaf 4 -beagle $subgenome.tmp.gz -out $pop.def.$subgenome -fai $HREF.fai -anc $HREF -sites anc.hg19.$subgenome.bed
		zcat $pop.beagle.corrected.gz | sed "s/chr//g" | gzip > $subgenome.tmp.gz
		angsd -doMaf 4 -beagle $subgenome.tmp.gz -out $pop.cor.$subgenome -fai $HREF.fai -anc $HREF -sites anc.hg19.$subgenome.bed
	done
	gzip -d -f *.$subgenome.mafs.gz
done

#for prop in 0.1 0.25 0.5 1.0
# do 
# zcat FIN.beagle.$prop.gz | sed "s/chr//g" | gzip > FIN.beagle.tmp.gz
# zcat YRI.beagle.$prop.gz | sed "s/chr//g" | gzip > YRI.beagle.tmp.gz
#  zcat FIN.beagle.corrected.$prop.gz | sed "s/chr//g" | gzip > FIN.beagle.corrected.tmp.gz
# zcat YRI.beagle.corrected.$prop.gz | sed "s/chr//g" | gzip > YRI.beagle.corrected.tmp.gz
##angsd -doMaf 4 -beagle FIN.beagle.tmp.gz -out FIN.def.$prop.$subgenome -fai $HREF.fai -anc $HREF -sites anc.hg19.$subgenome.bed
#angsd -doMaf 4 -beagle YRI.beagle.tmp.gz -out YRI.def.$prop.$subgenome -fai $HREF.fai -anc $HREF -sites anc.hg19.$subgenome.bed
#angsd -doMaf 4 -beagle FIN.beagle.corrected.tmp.gz -out FIN.cor.$prop.$subgenome -fai $HREF.fai -anc $HREF -sites anc.hg19.$subgenome.bed
#angsd -doMaf 4 -beagle YRI.beagle.corrected.tmp.gz -out YRI.cor.$prop.$subgenome -fai $HREF.fai -anc $HREF -sites anc.hg19.$subgenome.bed
# done

