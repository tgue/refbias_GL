#!/bin/bash

module load bioinfo-tools
module load samtools
module load bwa/0.7.17

INBAM=$1
CORES=$2

fbase=$(basename $INBAM)



samtools bam2fq $INBAM | awk '{switch (NR % 4) {case 1: print $1"X"; break; case 2: print substr($1, 1, length($1)/2); break; case 3: print "+"; break; case 0: print substr($1, 1, length($1)/2)} }'  > ${TMPDIR}/1.fastq

samtools bam2fq $INBAM | awk '{switch (NR % 4) {case 1: print $1"Y"; break; case 2: print substr($1, 1+length($1)/2, length($1)); break; case 3: print "+"; break; case 0: print substr($1, 1+length($1)/2, length($1))} }'  > ${TMPDIR}/2.fastq

cat ${TMPDIR}/1.fastq ${TMPDIR}/2.fastq > ${TMPDIR}/allsplit.fastq



bwa aln -l 16500 -n 0.01 -o 2 -t $CORES /proj/snic2020-2-10/1000AncientGenomes/ref_seqs/hs37d5.fa ${TMPDIR}/allsplit.fastq \
| bwa samse /proj/snic2020-2-10/1000AncientGenomes/ref_seqs/hs37d5.fa - ${TMPDIR}/allsplit.fastq  \
| samtools view -F 4 -h - \
| samtools sort -o /proj/snic2020-2-10/sllstore2017087/nobackup/private/reference_bias/refbias_GL/KGP/splitread_bams/$fbase -

samtools index /proj/snic2020-2-10/sllstore2017087/nobackup/private/reference_bias/refbias_GL/KGP/splitread_bams/$fbase
