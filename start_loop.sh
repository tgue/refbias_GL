#!/bin/bash

PROJ=snic2018-8-267 #snic2018-8-267   #p2018003 #snic2017-7-259
CORES=5
let CORES_ALLOC=CORES-1 ;
DEPTH=1.0
i=2.0X
p=0.2
div=50000

rm -rf fastngsadmix_results
mkdir fastngsadmix_results
rm admixture_results/*
rm qpadm_results/*
rm ngsadmix_results/*
rm pcangsd_results/*
#rm pcangsd_selection/*

for i in {1..50}
do
for p in 0.1 0.3 0.5 0.7 0.9
do
for div in 5000 10000 25000 50000
do
sbatch -A $PROJ -n $CORES_ALLOC -M snowy -p core -t 10-0:00:00  -J msp_${p}_${i}_${div} msp_gg_loop_fngsa_correctedGL.sh $p ${i} $div $DEPTH $CORES
done
done
done


