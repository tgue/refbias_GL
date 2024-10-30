#!/bin/bash

module load bioinfo-tools
module load ANGSD/0.933

cut -f 1,4 -d ' ' chip.selected.MAF0.2.nomiss.tped > sites.txt

sleep 61
angsd sites index sites.txt

ls /proj/snic2020-2-10/sllstore2017087/nobackup/private/reference_bias/refbias_GL/KGP/splitread_bams/*bam > tmp.bamfiles.txt

angsd -bam tmp.bamfiles.txt -checkBamHeaders 0 -nThreads 4 -doHaploCall 1 -doCounts 1 -doGeno -4 -doPost 2 -doPlink 2 -minMapQ 30 -minQ 30 -sites sites.txt -doMajorMinor 1 -GL 1 -domaf 1 -out pseudohaploid.all
haploToPlink pseudohaploid.all.haplo.gz pseudohaploid.all
sed -i "s/N/0/g" pseudohaploid.all.tped
