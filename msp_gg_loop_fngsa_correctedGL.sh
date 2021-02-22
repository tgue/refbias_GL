#!/bin/bash

module load bioinfo-tools
module load samtools/1.5_debug
module load bwa/0.7.17
module load AdapterRemoval/2.1.7
module load R/3.3.2


ITERATION=$2
SEQ_LEN=100000
SAMPLE=20
P1=0.0
P2=$1
NCHR=20
CORES=1
DIVERGENCE=$3
SEQ_COVERAGE=$4
PROJ_PATH=/proj/snic2020-2-10/sllstore2017087/nobackup/private/reference_bias/refbias_GL/
REFDIR=ref_seqs/
OUTDIR=sim_bams/

cp *sh *py $TMPDIR
cp -r cont/ $TMPDIR
cp -r endo/ $TMPDIR
cp -r gargammel/ $TMPDIR
cp -r data/ $TMPDIR

cd $TMPDIR

#create temp output dirs
mkdir ref_seqs
mkdir sim_fastq
mkdir sim_bam

python run_msprime_prep_gargammel.py --bases $SEQ_LEN --sample_per_pop $SAMPLE --seed $RANDOM --adm_prop $P1 $P2 --nchr $NCHR --divergence $DIVERGENCE

for refpop in S1 S2 S3
do
	cat endo/ancient.0.$refpop.1.*.fa > ref_seqs/ancient.$refpop.refgenome.fa
	bwa index ref_seqs/ancient.$refpop.refgenome.fa  
done

wait

for pop in T S1 S2 S3 O1 O2 O3 O4
	do
	for (( ind=1; ind<=$SAMPLE; ind++ ))
		do
		for (( c=1; c<=$NCHR; c++ ))
			do
			rm data/endo/*
			cp endo/ancient.$ind.$pop.1.chr${c}.fa data/endo/endo.1.fa
			cp endo/ancient.$ind.$pop.2.chr${c}.fa data/endo/endo.2.fa
			perl gargammel/gargammel.pl -c $SEQ_COVERAGE --comp 0,0,1 -l 40 -mapdamage gargammel/examplesMapDamage/results_LaBrana/misincorporation.txt double -o sim_fastq/ancient.$ind.$pop.chr${c} data/ 
			done

		zcat sim_fastq/ancient.$ind.$pop.chr*_s1.fq* | gzip > sim_fastq/ancient.$ind.$pop.allchr_s1.fq.gz 
		zcat sim_fastq/ancient.$ind.$pop.chr*_s2.fq* | gzip > sim_fastq/ancient.$ind.$pop.allchr_s2.fq.gz 

		filebase=ancient.$ind.$pop.allchr

		AdapterRemoval --file1 sim_fastq/ancient.$ind.$pop.allchr_s1.fq.gz --file2 sim_fastq/ancient.$ind.$pop.allchr_s2.fq.gz  --qualitybase 33 --gzip --qualitymax 60 --trimns --collapse --minalignmentlength 11 --threads $CORES --basename $TMPDIR/$filebase.merged --settings $TMPDIR/$filebase.out.settings

		cat $TMPDIR/$filebase.merged.collapsed.gz $TMPDIR/$filebase.merged.collapsed.truncated.gz $TMPDIR/$filebase.merged.pair1.truncated.gz $TMPDIR/$filebase.merged.pair2.truncated.gz > $TMPDIR/$filebase.merged.all.fastq.gz

		for REF in $(ls $REFDIR/*.fa)
			do
			only_ref=$(basename $REF)
			bwa aln -l 16500 -n 0.01 -o 2 -t $CORES $REF $TMPDIR/$filebase.merged.all.fastq.gz \
			| bwa samse $REF - $TMPDIR/$filebase.merged.all.fastq.gz  \
			| samtools view -F 4 -h -Su - \
			| samtools sort -o $TMPDIR/$filebase.$only_ref.bam -	
			samtools index $TMPDIR/$filebase.$only_ref.bam	
			#cp $TMPDIR/$filebase.$only_ref.bam.bai $TMPDIR/$filebase.$only_ref.bam ${PROJ_PATH}/BAMs
			done
		done
	done

wait


#for REF in $(ls ${REFDIR}/*.fa)
#	do
#	only_ref=$(basename $REF)
#	for file in $(ls *.${only_ref}.merged.bam)
#		do
#		filebase=$(basename ${file} .${only_ref}.merged.bam)
#		samtools view -F 4 -h ${file} | FilterUniqueSAMCons_rand.py | samtools view -h -Su - > ${filebase}.${only_ref}.merged.cons.bam 
#		done
#	wait
#	done

#wait


#for REF in $(ls ${REFDIR}/*.fa)
#	do
#	only_ref=$(basename $REF)
#	for file in $(ls *.${only_ref}.merged.bam)
#		do
#		filebase=$(basename ${file} .${only_ref}.merged.bam)
#		samtools calmd $TMPDIR/$filebase.$only_ref.merged.bam $REF | percidentity_threshold.py 0.9 35 $TMPDIR/$filebase.$only_ref.short.txt | samtools view -bS - > $filebase.$only_ref.cons.90perc.bam 
#		done
#	wait
#	done

#wait





#for file in $(ls *cons.90perc.bam)
#	do
#	samtools index $file 
#	done

#wait

for refpop in S1 S2 S3
	do
	#./beagleGL_fastngsadm.sh $refpop $P2 $ITERATION
	./beagleGL_fastngsadm_ngsadmref_correctedGL.sh $refpop $P2 $ITERATION $DIVERGENCE
	done

wait


