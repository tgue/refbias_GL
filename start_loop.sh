#!/bin/bash

PROJ=snic2017-7-259 #p2018003 #snic2017-7-259
CORES=1
i=1
p=0.71
div=10000

rm -rf fastngsadmix_results
mkdir fastngsadmix_results
rm admixture_results/*
rm qpadm_results/*

for i in {1..50}
do
for p in 0.1 0.3 0.5 0.7 0.9
do
for div in 5000 10000 25000 50000
do
sbatch -A $PROJ -n $CORES -p core -t 6-0:00:00  -J msp_${p}_${i}_${div} msp_gg_loop_fngsa_correctedGL.sh $p ${i} $div 1
done
done
done


