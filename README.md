# Estimating ancestry proportions and genotype likelihoods in the presence of mapping bias

This repository contains scripts used to conduct simulations of realistic aDNA sequence data from several populations under a demographic model. The data is then used as input for different ancestry estimation tools. Furthermore, we provide scripts to calculate a modified genoype likelihood for ascertained biallelic SNPs accounting for mapping bias.

(These scripts are far from being user-friendly or very efficient/optimized.)

## Calculating corrected genotype likelihoods

### Measuring mapping bias per individual and site
First, we need to assess mapping bias per individual. This step follows what we presented in [Günther and Nettelblad (2019)](https://doi.org/10.1371/journal.pgen.1008302) and (https://bitbucket.org/tguenther/refbias/src/master/).

Assuming that we have a Plink TPED file with biallelic SNPs and a reference genome, we prepare a list of SNPs with reference and alternative allele (using BioPython) and a corresponding .bed file:

```
python prepare_snplist.py $TPED $REFGENOME > snplist.txt

awk '{print $1,$2-1,$2}' snplist.txt > snplist.txt.bed
```

Next, we extract reads mapping to the SNP positions into a temporary BAM file and then modify each read at the SNP position, saving the results in a new BAM file:

```
samtools view -L snplist.txt.bed -h -b -u $BAMFILE > $BAMFILE.bedfilter.bam

samtools index	$BAMFILE.bedfilter.bam

samtools calmd $BAMFILE.bedfilter.bam $REFGENOME | ./modify_read_alternative.py snplist.txt | samtools view -b -u - > $BAMFILE.modified.bam
```

Now we remap those modified reads to the reference genome using the same mapping parameters that were used to produce the initial BAM file:

```
bwa aln -l 16500 -n 0.01 -o 2 -t $CORES $REFGENOME -b $BAMFILE.modified.bam \
| bwa samse $REFGENOME - $BAMFILE.modified.bam  \
| samtools view -F 4 -h - \
| samtools sort -o $BAMFILE.modified.remap.bam -

samtools index $BAMFILE.modified.remap.bam

samtools calmd -b $BAMFILE.modified.remap.bam $REFGENOME > $BAMFILE.modified.remap.MD.bam
samtools index $BAMFILE.modified.remap.MD.bam
```

Finally, we merge the new BAM file containing the modified reads with the bam file containing the original reads. the script `count_reads_mpB.py` is then counting the reads and outputs reference and alternative reads into `stdout`. Both numbers should be identical without mapping bias.

```
samtools merge $BAMFILE.merged.bam $BAMFILE.modified.remap.MD.bam $BAMFILE.bedfilter.bam


python count_reads_mpB.py $BAMFILE.merged.bam $REFGENOME $TPED snplist.txt.bed > counts.$BAMFILE
```

### GL calculation

Now we can calculate the actualy genotype likelihood using the python script `GL_calculator.py`. It requires the **original** BAM file (i.e. no modified reads), the list of SNPs, the counts file containing ref and alt allele counts, the reference genome and a name for the output file.
(Originally, we implemented this step to use the `pysam` libary but it turns out that parsing the `samtools mpileup` output was much faster, so that is the default version now.)

```
python GL_calculator.py $BAMFILE snplist.txt counts.$BAMFILE $REFGENOME output.gl
```

The script will produce two Beagle GL format output files `output.gl` and `output.gl.corrected` containing the default and corrected genotype likelihoods, respectively. If the files exist already, additional columns are added to the files. By looping through a number of individuals and their BAM and count files, one can create a GL file for multiple individuals/populations.

## Simulations

The folder `simulations` contains scripts for running the simulations for the study, using `msprime` and `gargammel` to generate the data and then `ADMIXTURE`, `qpAdm`, `ngsadmix` and `fastngsadmix` to estimate ancestry proportions. The simulations are started from the script `msp_gg_loop_fngsa_correctedGL.sh`.
