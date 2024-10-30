#!/bin/bash

module load bioinfo-tools
module load samtools
module load bwa/0.7.17

INBAM=$1
CORES=$2

fbase=$(basename $INBAM)

samtools view -L snps.bed -h -b -u -q 30 $INBAM > ${TMPDIR}/bedfilter.bam

samtools index ${TMPDIR}/bedfilter.bam
samtools view -c ${TMPDIR}/bedfilter.bam

#samtools bam2fq ${TMPDIR}/bedfilter.bam > ${TMPDIR}/bedfilter.fastq

samtools collate -u -O ${TMPDIR}/bedfilter.bam | samtools fastq -1 ${TMPDIR}/$c.paired1.fq -2 ${TMPDIR}/$c.paired2.fq -0 /dev/null -s ${TMPDIR}/$c.single.fq -n

cat ${TMPDIR}/$c.paired1.fq ${TMPDIR}/$c.paired2.fq ${TMPDIR}/$c.single.fq > ${TMPDIR}/bedfilter.fastq

wc -l ${TMPDIR}/bedfilter.fastq

bwa aln -l 16500 -n 0.01 -o 2 -t $CORES /proj/snic2020-2-10/1000AncientGenomes/ref_seqs/hs37d5.fa ${TMPDIR}/bedfilter.fastq \
| bwa samse /proj/snic2020-2-10/1000AncientGenomes/ref_seqs/hs37d5.fa - ${TMPDIR}/bedfilter.fastq  \
| samtools view -F 4 -h - \
| samtools sort -o /proj/snic2020-2-10/sllstore2017087/nobackup/private/reference_bias/refbias_GL/KGP/remapped_bams/$fbase -

samtools index /proj/snic2020-2-10/sllstore2017087/nobackup/private/reference_bias/refbias_GL/KGP/remapped_bams/$fbase
