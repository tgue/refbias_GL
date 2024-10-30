#!/bin/bash

module load bioinfo-tools plink/1.90b4.9

cut -d '/' -f 7 bamlist.txt > selected_IDs.txt

wget  ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.vcf.gz

plink --vcf ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.vcf.gz --snps-only --make-bed --allow-extra-chr --out  ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes

grep -f selected_IDs.txt ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.fam > keep.txt

plink --bfile ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes --snps-only --keep keep.txt --make-bed --out chip.selected --autosome --allow-extra-chr

plink --bfile chip.selected --keep keep.txt --maf 0.2 --recode transpose --geno 0 --out chip.selected.MAF0.2.nomiss

#manually edit chip.selected.MAF0.2.nomiss.tfam to add populations

python tped2beagleGL.py FIN
python tped2beagleGL.py YRI
python tped2beagleGL.py JPT

gzip FIN* YRI* JPT*

rm ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.*
