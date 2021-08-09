#!/bin/bash

module load bioinfo-tools
module load ANGSD/0.917
module load plink/1.90b4.9
module load R/3.3.2
module load samtools/1.5_debug
module load bwa/0.7.17
module load ADMIXTURE/1.3.0
module load AdmixTools/7.0.1
module load python/3.8.7
module load pysam
module load biopython

REFPOP=$1
PROP=$2
ITERATION=$3
DIVERGENCE=$4
CORES=$5
REFERENCE_SAMPLESIZE=20
TARGET_SAMPLESIZE=20
PROJ_PATH=/proj/snic2020-2-10/sllstore2017087/nobackup/private/reference_bias/refbias_GL/

ls *.S1.allchr.ancient.S1.refgenome.fa.bam > $REFPOP.$ITERATION.bamlist

#ascertain SNPs in outgroup
angsd -GL 1 -out ${TMPDIR}/$REFPOP.$ITERATION.genolike -nThreads $CORES -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-7 -doMaf 1 -minMaf 0.1 -bam $REFPOP.$ITERATION.bamlist -minQ 30 -minMapQ 30 -remove_bads 1 -uniqueOnly 1 -rmTrans 1

#select maximum 200000 SNPs and format SNP list
zcat ${TMPDIR}/$REFPOP.$ITERATION.genolike.beagle.gz | tail -n +2 | cut -f 1 | shuf -n 200000 >  ${TMPDIR}/$REFPOP.sites.$ITERATION.tmp
cut -f 1 -d '_' ${TMPDIR}/$REFPOP.sites.$ITERATION.tmp > ${TMPDIR}/$REFPOP.chr.$ITERATION.tmp
cut -f 2 -d '_' ${TMPDIR}/$REFPOP.sites.$ITERATION.tmp > ${TMPDIR}/$REFPOP.bp.$ITERATION.tmp
paste ${TMPDIR}/$REFPOP.chr.$ITERATION.tmp  ${TMPDIR}/$REFPOP.bp.$ITERATION.tmp  | sort -b -k 1,1 -k 2,2n  > ${TMPDIR}/$REFPOP.$ITERATION.sites.txt
#cp ${TMPDIR}/$REFPOP.$ITERATION.genolike.beagle.gz $PROJ_PATH
#cp ${TMPDIR}/$REFPOP.$ITERATION.sites.txt $PROJ_PATH

ls *.T.allchr.ancient.$REFPOP.refgenome.fa.bam > ${TMPDIR}/$REFPOP.$ITERATION.T.bamlist
rm ${TMPDIR}/$REFPOP.$ITERATION.sites.txt.idx ${TMPDIR}/$REFPOP.$ITERATION.sites.txt.idx

sleep 61
angsd sites index ${TMPDIR}/$REFPOP.$ITERATION.sites.txt

##prepare for pruning

#do haplocalls in all
ls *.T.allchr.ancient.$REFPOP.refgenome.fa.bam *.S2.allchr.ancient.$REFPOP.refgenome.fa.bam *.S3.allchr.ancient.$REFPOP.refgenome.fa.bam | rev | sort | rev > ${TMPDIR}/$REFPOP.$ITERATION.all.bamlist
cp ${TMPDIR}/$REFPOP.$ITERATION.all.bamlist ${TMPDIR}/$REFPOP.$ITERATION.noO.bamlist
angsd -bam ${TMPDIR}/$REFPOP.$ITERATION.all.bamlist -checkBamHeaders 0 -nThreads $CORES -doHaploCall 1 -doCounts 1 -doGeno -4 -doPost 2 -doPlink 2 -minMapQ 30 -minQ 30 -sites ${TMPDIR}/$REFPOP.$ITERATION.sites.txt -doMajorMinor 1 -GL 1 -domaf 1 -out ${TMPDIR}/$REFPOP.$ITERATION.all


ls *.S2.allchr.ancient.$REFPOP.refgenome.fa.bam *.S3.allchr.ancient.$REFPOP.refgenome.fa.bam | rev | sort | rev > ${TMPDIR}/$REFPOP.$ITERATION.sources.bamlist

#cp ${TMPDIR}/$REFPOP.$ITERATION.all.bamlist $PROJ_PATH

plink --tfile ${TMPDIR}/$REFPOP.$ITERATION.all --list-duplicate-vars --allow-extra-chr
plink --tfile ${TMPDIR}/$REFPOP.$ITERATION.all  --geno 0.5 --exclude plink.dupvar --alleleACGT --allow-extra-chr --make-bed --out ${TMPDIR}/$REFPOP.$ITERATION.all
plink --tfile ${TMPDIR}/$REFPOP.$ITERATION.all  --geno 0.5 --exclude plink.dupvar --allow-extra-chr --alleleACGT --recode transpose --out ${TMPDIR}/$REFPOP.$ITERATION.all

cut -f 1,4,5,6 ${TMPDIR}/$REFPOP.$ITERATION.all.bim > ${TMPDIR}/$REFPOP.$ITERATION.sites.txt
sed -i "s/^/chr/g" ${TMPDIR}/$REFPOP.$ITERATION.sites.txt
sed -i "s/^/chr/g" ${TMPDIR}/$REFPOP.$ITERATION.all.tped

python prepare_snplist.py ${TMPDIR}/$REFPOP.$ITERATION.all.tped ref_seqs/ancient.$REFPOP.refgenome.fa > $REFPOP.refgenome.fa.snplist.nonpruned.txt

##pruning

plink --bfile ${TMPDIR}/$REFPOP.$ITERATION.all --indep-pairwise 200 25 0.4
plink --bfile ${TMPDIR}/$REFPOP.$ITERATION.all --extract plink.prune.in --make-bed --out ${TMPDIR}/$REFPOP.$ITERATION.all.pruned
plink --bfile ${TMPDIR}/$REFPOP.$ITERATION.all.pruned --recode transpose --out ${TMPDIR}/$REFPOP.$ITERATION.all.pruned

cut -f 1,4,5,6 ${TMPDIR}/$REFPOP.$ITERATION.all.pruned.bim > ${TMPDIR}/$REFPOP.$ITERATION.sites.pruned.txt
sed -i "s/^/chr/g" ${TMPDIR}/$REFPOP.$ITERATION.sites.pruned.txt
sed -i "s/^/chr/g" ${TMPDIR}/$REFPOP.$ITERATION.all.pruned.tped
#cp ${TMPDIR}/$REFPOP.$ITERATION.sites.pruned.txt $PROJ_PATH

sleep 61

angsd sites index ${TMPDIR}/$REFPOP.$ITERATION.sites.txt
angsd sites index ${TMPDIR}/$REFPOP.$ITERATION.sites.pruned.txt

sleep 61

#do haplocalls in all including outgroups for qpadm
ls *.O*.allchr.ancient.$REFPOP.refgenome.fa.bam | rev | sort | rev >> ${TMPDIR}/$REFPOP.$ITERATION.all.bamlist
angsd -bam ${TMPDIR}/$REFPOP.$ITERATION.all.bamlist -checkBamHeaders 0 -nThreads $CORES -doHaploCall 1 -doCounts 1 -doGeno -4 -doPost 2 -doPlink 2 -minMapQ 30 -minQ 30 -sites ${TMPDIR}/$REFPOP.$ITERATION.sites.txt -doMajorMinor 1 -GL 1 -domaf 1 -out ${TMPDIR}/$REFPOP.$ITERATION.allpop

rm *.O*.allchr.ancient.$REFPOP.refgenome.fa.bam
#cp ${TMPDIR}/$REFPOP.$ITERATION.all.bamlist $PROJ_PATH

plink --tfile ${TMPDIR}/$REFPOP.$ITERATION.allpop --list-duplicate-vars --allow-extra-chr
plink --tfile ${TMPDIR}/$REFPOP.$ITERATION.allpop --geno 0.5 --exclude plink.dupvar --allow-extra-chr --alleleACGT --make-bed --out ${TMPDIR}/allpop

plink --bfile ${TMPDIR}/allpop --extract plink.prune.in --allow-extra-chr --alleleACGT --make-bed --out ${TMPDIR}/allpop.pruned

#cp ${TMPDIR}/$REFPOP.$ITERATION.all* $PROJ_PATH



#prep for modified reads analysis

python prepare_snplist.py ${TMPDIR}/$REFPOP.$ITERATION.all.tped ref_seqs/ancient.$REFPOP.refgenome.fa > $REFPOP.refgenome.fa.snplist.txt


#modify reads and count in all
index=1
cat ${TMPDIR}/$REFPOP.$ITERATION.noO.bamlist | while read bamfile; do
	#modify reads
	samtools calmd $bamfile ref_seqs/ancient.$REFPOP.refgenome.fa | ./modify_read_alternative.py $REFPOP.refgenome.fa.snplist.txt | samtools view -b -u - > modified.$bamfile

	#remap
	bwa aln -l 16500 -n 0.01 -o 2 -t $CORES ref_seqs/ancient.$REFPOP.refgenome.fa -b modified.$bamfile \
	| bwa samse ref_seqs/ancient.$REFPOP.refgenome.fa - modified.$bamfile  \
	| samtools view -F 4 -h - \
	| samtools sort -o modified.remap.$bamfile -

	samtools index modified.remap.$bamfile

	samtools calmd -b modified.remap.$bamfile ref_seqs/ancient.$REFPOP.refgenome.fa > modified.remap.MD.$bamfile

	samtools index modified.remap.MD.$bamfile

	samtools merge merged.modified.remap.$bamfile modified.remap.MD.$bamfile $bamfile

	awk '{print $1,$2-1,$2}' $REFPOP.refgenome.fa.snplist.txt > $REFPOP.refgenome.fa.snplist.txt.bed
	python count_reads_mpB.py merged.modified.remap.$bamfile ref_seqs/ancient.$REFPOP.refgenome.fa ${TMPDIR}/$REFPOP.$ITERATION.all.pruned.tped $REFPOP.refgenome.fa.snplist.txt.bed > counts.$bamfile
done


angsd sites index ${TMPDIR}/$REFPOP.$ITERATION.sites.pruned.txt

sleep 10

###this step is using pruned sites###

#calculate GL in T
cat ${TMPDIR}/$REFPOP.$ITERATION.T.bamlist | while read bamfile; do
python pysam_worker.py $bamfile ${TMPDIR}/$REFPOP.$ITERATION.sites.txt counts.$bamfile ref_seqs/ancient.$REFPOP.refgenome.fa ${TMPDIR}/$REFPOP.$ITERATION.T.gl
done


#calculate GL in sources
cat ${TMPDIR}/$REFPOP.$ITERATION.sources.bamlist | while read bamfile; do
python pysam_worker.py $bamfile ${TMPDIR}/$REFPOP.$ITERATION.sites.txt counts.$bamfile ref_seqs/ancient.$REFPOP.refgenome.fa ${TMPDIR}/$REFPOP.$ITERATION.sources.beagle
done
cp ${TMPDIR}/$REFPOP.$ITERATION.sources.beagle ${TMPDIR}/$REFPOP.$ITERATION.sources.beagle.input
cp ${TMPDIR}/$REFPOP.$ITERATION.sources.beagle.corrected ${TMPDIR}/$REFPOP.$ITERATION.sources.beagle.input.corrected


#preparation for fastNGSadmix
#run NGSadmix on uncorrected sources

NGSadmix -likes ${TMPDIR}/$REFPOP.$ITERATION.sources.beagle.input -P $CORES -K 2 -printInfo 1 -outfiles ${TMPDIR}/$REFPOP.$ITERATION.sources.NGSADM.beagle

#cp ${TMPDIR}/$REFPOP.$ITERATION.sources.NGSADM.beagle.* ${PROJ_PATH}
gunzip ${TMPDIR}/$REFPOP.$ITERATION.sources.NGSADM.beagle.fopt.gz

q1=$(head -n1 ${TMPDIR}/$REFPOP.$ITERATION.sources.NGSADM.beagle.qopt | cut -f1 -d ' ')

if [ $(echo $q1 '<' 0.4 | bc -l) -eq 1 ]
then
	awk '{print $2,$1}' ${TMPDIR}/$REFPOP.$ITERATION.sources.NGSADM.beagle.fopt > tmp.fopt
	cp tmp.fopt ${TMPDIR}/$REFPOP.$ITERATION.sources.NGSADM.beagle.fopt
fi


echo "id chr pos name A0_freq A1 S2 S3" > refPanel_${REFPOP}.$ITERATION.sources.txt

cut -f1,2,3 ${TMPDIR}/$REFPOP.$ITERATION.sources.beagle > tmp.beagle


Rscript ${PROJ_PATH}/makeRefNGSadmix.R tmp.beagle ${TMPDIR}/$REFPOP.$ITERATION.sources.NGSADM.beagle.filter
paste tmp.ref ${TMPDIR}/$REFPOP.$ITERATION.sources.NGSADM.beagle.fopt >> refPanel_${REFPOP}.$ITERATION.sources.txt
rm tmp.ref tmp.beagle

echo S2 S3 >> ${TMPDIR}/$REFPOP.$ITERATION.fngs_ref_ninds.txt
NREF_HAPL=40
echo $NREF_HAPL $NREF_HAPL >> ${TMPDIR}/$REFPOP.$ITERATION.fngs_ref_ninds.txt


#run NGSadmix on corrected sources

NGSadmix -likes ${TMPDIR}/$REFPOP.$ITERATION.sources.beagle.input.corrected -P $CORES -K 2 -printInfo 1 -outfiles ${TMPDIR}/$REFPOP.$ITERATION.sources.NGSADM.beagle.corrected

gunzip ${TMPDIR}/$REFPOP.$ITERATION.sources.NGSADM.beagle.corrected.fopt.gz

q1=$(head -n1 ${TMPDIR}/$REFPOP.$ITERATION.sources.NGSADM.beagle.corrected.qopt | cut -f1 -d ' ')

if [ $(echo $q1 '<' 0.4 | bc -l) -eq 1 ]
then
	awk '{print $2,$1}' ${TMPDIR}/$REFPOP.$ITERATION.sources.NGSADM.beagle.corrected.fopt > tmp.fopt
	cp tmp.fopt ${TMPDIR}/$REFPOP.$ITERATION.sources.NGSADM.beagle.corrected.fopt
fi


echo "id chr pos name A0_freq A1 S2 S3" > refPanel_${REFPOP}.$ITERATION.sources.corrected.txt

cut -f1,2,3 ${TMPDIR}/$REFPOP.$ITERATION.sources.beagle.input.corrected > tmp.beagle


Rscript ${PROJ_PATH}/makeRefNGSadmix.R tmp.beagle ${TMPDIR}/$REFPOP.$ITERATION.sources.NGSADM.beagle.corrected.filter
paste tmp.ref ${TMPDIR}/$REFPOP.$ITERATION.sources.NGSADM.beagle.corrected.fopt >> refPanel_${REFPOP}.$ITERATION.sources.corrected.txt
rm tmp.ref tmp.beagle


#run fastNGSadmix per individual

#prune GL files
head -n1 ${TMPDIR}/$REFPOP.$ITERATION.T.gl > gl.head
cat ${TMPDIR}/$REFPOP.$ITERATION.T.gl | python2 prune_gl.py plink.prune.in > tmp.gl
cat gl.head tmp.gl > ${TMPDIR}/$REFPOP.$ITERATION.T.pruned.gl

cat ${TMPDIR}/$REFPOP.$ITERATION.T.gl.corrected | python2 prune_gl.py plink.prune.in > tmp.gl
cat gl.head tmp.gl > ${TMPDIR}/$REFPOP.$ITERATION.T.pruned.gl.corrected

index=1
cat ${TMPDIR}/$REFPOP.$ITERATION.T.bamlist | while read line; do
python ${PROJ_PATH}/uniq_gl.py ${TMPDIR}/$REFPOP.$ITERATION.T.gl $index refPanel_${REFPOP}.$ITERATION.sources.txt
python ${PROJ_PATH}/uniq_gl.py ${TMPDIR}/$REFPOP.$ITERATION.T.gl.corrected $index refPanel_${REFPOP}.$ITERATION.sources.txt
python ${PROJ_PATH}/uniq_gl.py ${TMPDIR}/$REFPOP.$ITERATION.T.pruned.gl $index refPanel_${REFPOP}.$ITERATION.sources.txt
python ${PROJ_PATH}/uniq_gl.py ${TMPDIR}/$REFPOP.$ITERATION.T.pruned.gl.corrected $index refPanel_${REFPOP}.$ITERATION.sources.txt

fastNGSadmix -likes $REFPOP.$ITERATION.T.gl.uniq -fname refPanel_${REFPOP}.$ITERATION.sources.txt -Nname ${TMPDIR}/$REFPOP.$ITERATION.fngs_ref_ninds.txt -out /proj/snic2020-2-10/sllstore2017087/nobackup/private/reference_bias/refbias_GL/fastngsadmix_results/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.ngsadmix_beagleGL_ref${REFPOP}_${index} -whichPops S2,S3 
fastNGSadmix -likes $REFPOP.$ITERATION.T.gl.corrected.uniq -fname refPanel_${REFPOP}.$ITERATION.sources.corrected.txt -Nname ${TMPDIR}/$REFPOP.$ITERATION.fngs_ref_ninds.txt -out /proj/snic2020-2-10/sllstore2017087/nobackup/private/reference_bias/refbias_GL/fastngsadmix_results/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.ngsadmix_beagleGL_ref${REFPOP}_${index}.corrected -whichPops S2,S3 
fastNGSadmix -likes $REFPOP.$ITERATION.T.pruned.gl.uniq -fname refPanel_${REFPOP}.$ITERATION.sources.txt -Nname ${TMPDIR}/$REFPOP.$ITERATION.fngs_ref_ninds.txt -out /proj/snic2020-2-10/sllstore2017087/nobackup/private/reference_bias/refbias_GL/fastngsadmix_results/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.ngsadmix_beagleGL_ref${REFPOP}_${index}.pruned -whichPops S2,S3 
fastNGSadmix -likes $REFPOP.$ITERATION.T.pruned.gl.corrected.uniq -fname refPanel_${REFPOP}.$ITERATION.sources.corrected.txt -Nname ${TMPDIR}/$REFPOP.$ITERATION.fngs_ref_ninds.txt -out /proj/snic2020-2-10/sllstore2017087/nobackup/private/reference_bias/refbias_GL/fastngsadmix_results/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.ngsadmix_beagleGL_ref${REFPOP}_${index}.corrected.pruned -whichPops S2,S3 

let "index += 1" 
done


yes "S2" | head -n $REFERENCE_SAMPLESIZE > ${TMPDIR}/$REFPOP.$ITERATION.all.pruned.pop
yes "S3" | head -n $REFERENCE_SAMPLESIZE >> ${TMPDIR}/$REFPOP.$ITERATION.all.pruned.pop
yes "-" | head -n $TARGET_SAMPLESIZE >> ${TMPDIR}/$REFPOP.$ITERATION.all.pruned.pop

cp ${TMPDIR}/$REFPOP.$ITERATION.all.pruned.pop ${TMPDIR}/$REFPOP.$ITERATION.all.pop


admixture --haploid='*' --supervised -s $RANDOM ${TMPDIR}/$REFPOP.$ITERATION.all.pruned.bed 2
mv ${TMPDIR}/$REFPOP.$ITERATION.all.pruned.2.Q /proj/snic2020-2-10/sllstore2017087/nobackup/private/reference_bias/refbias_GL/admixture_results/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.supervised.admixture.2.Q

admixture --haploid='*' --supervised -s $RANDOM ${TMPDIR}/$REFPOP.$ITERATION.all.bed 2
mv ${TMPDIR}/$REFPOP.$ITERATION.all.2.Q /proj/snic2020-2-10/sllstore2017087/nobackup/private/reference_bias/refbias_GL/admixture_results/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.supervised.admixture.nonpruned.2.Q

admixture --haploid='*' --supervised -j${CORES} -s $RANDOM ${TMPDIR}/$REFPOP.$ITERATION.all.pruned.bed 2
mv ${TMPDIR}/$REFPOP.$ITERATION.all.pruned.2.Q /proj/snic2020-2-10/sllstore2017087/nobackup/private/reference_bias/refbias_GL/admixture_results/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.supervised.multithread.admixture.2.Q

admixture --haploid='*' --supervised -j${CORES} -s $RANDOM ${TMPDIR}/$REFPOP.$ITERATION.all.bed 2
mv ${TMPDIR}/$REFPOP.$ITERATION.all.2.Q /proj/snic2020-2-10/sllstore2017087/nobackup/private/reference_bias/refbias_GL/admixture_results/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.supervised.multithread.admixture.nonpruned.2.Q

admixture --supervised -j${CORES} -s $RANDOM ${TMPDIR}/$REFPOP.$ITERATION.all.pruned.bed 2
mv ${TMPDIR}/$REFPOP.$ITERATION.all.pruned.2.Q /proj/snic2020-2-10/sllstore2017087/nobackup/private/reference_bias/refbias_GL/admixture_results/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.supervised.nonhaploid.admixture.2.Q

admixture --supervised -j${CORES} -s $RANDOM ${TMPDIR}/$REFPOP.$ITERATION.all.bed 2
mv ${TMPDIR}/$REFPOP.$ITERATION.all.2.Q /proj/snic2020-2-10/sllstore2017087/nobackup/private/reference_bias/refbias_GL/admixture_results/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.supervised.nonhaploid.admixture.nonpruned.2.Q

admixture --haploid='*' -s $RANDOM ${TMPDIR}/$REFPOP.$ITERATION.all.pruned.bed 2
mv ${TMPDIR}/$REFPOP.$ITERATION.all.pruned.2.Q /proj/snic2020-2-10/sllstore2017087/nobackup/private/reference_bias/refbias_GL/admixture_results/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.unsupervised.admixture.2.Q

admixture --haploid='*' -s $RANDOM ${TMPDIR}/$REFPOP.$ITERATION.all.bed 2
mv ${TMPDIR}/$REFPOP.$ITERATION.all.2.Q /proj/snic2020-2-10/sllstore2017087/nobackup/private/reference_bias/refbias_GL/admixture_results/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.unsupervised.admixture.nonpruned.2.Q

##qpAdm

yes "S2" | head -n $REFERENCE_SAMPLESIZE > ${TMPDIR}/all.tmp.pop
yes "S3" | head -n $REFERENCE_SAMPLESIZE >> ${TMPDIR}/all.tmp.pop
yes "T" | head -n $TARGET_SAMPLESIZE >> ${TMPDIR}/all.tmp.pop
yes "O1" | head -n $REFERENCE_SAMPLESIZE >> ${TMPDIR}/all.tmp.pop
yes "O2" | head -n $REFERENCE_SAMPLESIZE >> ${TMPDIR}/all.tmp.pop
yes "O3" | head -n $REFERENCE_SAMPLESIZE >> ${TMPDIR}/all.tmp.pop
yes "O4" | head -n $REFERENCE_SAMPLESIZE >> ${TMPDIR}/all.tmp.pop

rm ${TMPDIR}/all.index.ind

NUM_INDS=$(($TARGET_SAMPLESIZE+6*$REFERENCE_SAMPLESIZE))

for (( i=1; i<=$NUM_INDS; i++ ))
do
   echo $i >> ${TMPDIR}/all.index.ind
done

cut -d ' ' -f 3,4,5 ${TMPDIR}/allpop.fam > ${TMPDIR}/ind.tmp

paste -d ' ' ${TMPDIR}/all.tmp.pop ${TMPDIR}/all.index.ind ${TMPDIR}/ind.tmp  ${TMPDIR}/all.tmp.pop > ${TMPDIR}/allpop.fam

cp ${PROJ_PATH}/*.par ${PROJ_PATH}/*.pops ${TMPDIR}

#non-pruned

convertf -p convertf.par

qpAdm -p qpadm.par > $PROJ_PATH/qpadm_results/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.qpadm.out

#pruned

mv allpop.pruned.bed allpop.bed
mv allpop.pruned.bim allpop.bim

convertf -p convertf.par

qpAdm -p qpadm.par > $PROJ_PATH/qpadm_results/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.pruned.qpadm.out


##NGSadmix

ls *.T.allchr.ancient.$REFPOP.refgenome.fa.bam *.S2.allchr.ancient.$REFPOP.refgenome.fa.bam *.S3.allchr.ancient.$REFPOP.refgenome.fa.bam | grep -v counts | grep -v modified | rev | sort | rev > ${TMPDIR}/$REFPOP.$ITERATION.all.bamlist

rm ${TMPDIR}/$REFPOP.$ITERATION.all.beagle*

#calculate GL in all
cat ${TMPDIR}/$REFPOP.$ITERATION.all.bamlist | while read bamfile; do
python pysam_worker.py $bamfile ${TMPDIR}/$REFPOP.$ITERATION.sites.txt counts.$bamfile ref_seqs/ancient.$REFPOP.refgenome.fa ${TMPDIR}/$REFPOP.$ITERATION.all.beagle
done


NGSadmix -likes ${TMPDIR}/$REFPOP.$ITERATION.all.beagle -P $CORES -K 2 -outfiles ${TMPDIR}/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.all.NGSADM.beagle
cp ${TMPDIR}/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.all.NGSADM.beagle.qopt ${PROJ_PATH}/ngsadmix_results/
rm ${TMPDIR}/$REFPOP.$ITERATION.all.beagle.filter.gz
gzip ${TMPDIR}/$REFPOP.$ITERATION.all.beagle


NGSadmix -likes ${TMPDIR}/$REFPOP.$ITERATION.all.beagle.corrected -P $CORES -K 2 -outfiles ${TMPDIR}/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.all.NGSADM.beagle.corrected
cp ${TMPDIR}/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.all.NGSADM.beagle.corrected.qopt ${PROJ_PATH}/ngsadmix_results/
cp ${TMPDIR}/$REFPOP.$ITERATION.all.beagle.corrected ${TMPDIR}/$REFPOP.$ITERATION.all.beagle.corr
rm ${TMPDIR}/$REFPOP.$ITERATION.all.beagle.corr.gz
gzip ${TMPDIR}/$REFPOP.$ITERATION.all.beagle.corrected



#calculate GL in all
cat ${TMPDIR}/$REFPOP.$ITERATION.all.bamlist | while read bamfile; do
python pysam_worker.py $bamfile ${TMPDIR}/$REFPOP.$ITERATION.sites.pruned.txt counts.$bamfile ref_seqs/ancient.$REFPOP.refgenome.fa ${TMPDIR}/$REFPOP.$ITERATION.all.beagle
done


NGSadmix -likes ${TMPDIR}/$REFPOP.$ITERATION.all.beagle -P $CORES -K 2 -outfiles ${TMPDIR}/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.all.NGSADM.beagle.pruned
cp ${TMPDIR}/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.all.NGSADM.beagle.pruned.qopt ${PROJ_PATH}/ngsadmix_results/
cp ${TMPDIR}/$REFPOP.$ITERATION.all.beagle ${TMPDIR}/$REFPOP.$ITERATION.all.beagle.pruned

gzip ${TMPDIR}/$REFPOP.$ITERATION.all.beagle.pruned

NGSadmix -likes ${TMPDIR}/$REFPOP.$ITERATION.all.beagle.corrected -P $CORES -K 2 -outfiles ${TMPDIR}/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.all.NGSADM.beagle.pruned.corrected
cp ${TMPDIR}/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.all.NGSADM.beagle.pruned.corrected.qopt ${PROJ_PATH}/ngsadmix_results/
cp ${TMPDIR}/$REFPOP.$ITERATION.all.beagle.corrected ${TMPDIR}/$REFPOP.$ITERATION.all.beagle.corr.pruned

gzip ${TMPDIR}/$REFPOP.$ITERATION.all.beagle.corr.pruned




##PCAngsd

python /home/torsteng/pcangsd/pcangsd.py -beagle ${TMPDIR}/$REFPOP.$ITERATION.all.beagle.gz -admix -admix_K 2 -threads $CORES
python printQ.py
cp out.Q.2.txt ${PROJ_PATH}/pcangsd_results/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.all.pcangsd.uncorr.qopt

python /home/torsteng/pcangsd/pcangsd.py -beagle ${TMPDIR}/$REFPOP.$ITERATION.all.beagle.corrected.gz -admix -admix_K 2 -threads $CORES
python printQ.py
cp out.Q.2.txt ${PROJ_PATH}/pcangsd_results/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.all.pcangsd.corr.qopt

python /home/torsteng/pcangsd/pcangsd.py -beagle ${TMPDIR}/$REFPOP.$ITERATION.all.beagle.pruned.gz -admix -admix_K 2 -threads $CORES
python printQ.py
cp out.Q.2.txt ${PROJ_PATH}/pcangsd_results/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.all.pcangsd.uncorr.pruned.qopt

python /home/torsteng/pcangsd/pcangsd.py -beagle ${TMPDIR}/$REFPOP.$ITERATION.all.beagle.corr.pruned.gz -admix -admix_K 2 -threads $CORES
python printQ.py
cp out.Q.2.txt ${PROJ_PATH}/pcangsd_results/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.all.pcangsd.corr.pruned.qopt


##PCAngsd selection (only on pruned data for now)
module load R/3.6.1
module load R_packages/3.6.1

#python /home/torsteng/pcangsd/pcangsd.py -beagle ${TMPDIR}/$REFPOP.$ITERATION.all.beagle.pruned.gz -pcadapt -threads $CORES
#Rscript /home/torsteng/pcangsd/scripts/pcadapt.R pcangsd.pcadapt.zscores.npy
#cp NA.pcadapt.pval.txt ${PROJ_PATH}/pcangsd_selection/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.all.pcadapt.uncorr.pval.txt
#cp NA.pcadapt.test.txt ${PROJ_PATH}/pcangsd_selection/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.all.pcadapt.uncorr.test.txt

#python /home/torsteng/pcangsd/pcangsd.py -beagle ${TMPDIR}/$REFPOP.$ITERATION.all.beagle.corr.pruned.gz -pcadapt -threads $CORES
#Rscript /home/torsteng/pcangsd/scripts/pcadapt.R pcangsd.pcadapt.zscores.npy
#cp NA.pcadapt.pval.txt ${PROJ_PATH}/pcangsd_selection/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.all.pcadapt.corr.pval.txt
#cp NA.pcadapt.test.txt ${PROJ_PATH}/pcangsd_selection/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.all.pcadapt.corr.test.txt


#python /home/torsteng/pcangsd/pcangsd.py -beagle ${TMPDIR}/$REFPOP.$ITERATION.all.beagle.pruned.gz -selection -threads $CORES -selection_e 5 -sites_save
#python printSel.py pcangsd.selection.npy
#cp out.selection.2.txt ${PROJ_PATH}/pcangsd_selection/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.all.selection.uncorr.test.txt
#cp pcangsd.sites ${PROJ_PATH}/pcangsd_selection/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.all.selection.uncorr.sites.txt

#python /home/torsteng/pcangsd/pcangsd.py -beagle ${TMPDIR}/$REFPOP.$ITERATION.all.beagle.corr.pruned.gz -selection -threads $CORES -selection_e 5 -sites_save
#python printSel.py pcangsd.selection.npy
#cp out.selection.2.txt ${PROJ_PATH}/pcangsd_selection/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.all.selection.corr.test.txt
#cp pcangsd.sites ${PROJ_PATH}/pcangsd_selection/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.all.selection.corr.sites.txt



#cp ${TMPDIR}/$REFPOP.$ITERATION.all.beagle.pruned.gz ${PROJ_PATH}/pcangsd_selection/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.beagle.pruned.gz
#cp ${TMPDIR}/$REFPOP.$ITERATION.all.beagle.corr.pruned.gz ${PROJ_PATH}/pcangsd_selection/$REFPOP.$PROP.$DIVERGENCE.$ITERATION.beagle.corr.pruned.gz
#cp counts.* ${PROJ_PATH}/count_files/


rm /proj/snic2020-2-10/sllstore2017087/nobackup/private/reference_bias/refbias_GL/fastngsadmix_results/*log
