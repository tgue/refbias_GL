#!/bin/bash

module load bioinfo-tools
module load ANGSD/0.933

for pop in FIN YRI JPT
do
angsd -doMaf 4 -beagle $pop.pseudohaploid.gl.gz -out $pop.ph -fai /proj/snic2020-2-10/1000AncientGenomes/ref_seqs/hs37d5.fa.fai -anc /proj/snic2020-2-10/1000AncientGenomes/ref_seqs/hs37d5.fa

gzip -d $pop.ph.mafs.gz
done


