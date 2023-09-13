#!/bin/bash
#SBATCH --job-name="plink"
#SBATCH --export=ALL
#SBATCH --mem=90G

# PLINK

# Follow tutorial:
# https://zzz.bwh.harvard.edu/plink/tutorial.shtml

#path="/home/tmichel/scratch/MINIMAP_long_reads/data" ;

mkdir ./plink

#species="coelocentrum"

# Estimate ROHi for gVCF file:
plink \
 --vcf GVCF_SNPs.vcf.gz \
 --homozyg \
 --out ./plink/GVCF_SNPs \
 --allow-extra-chr \
 --no-parents \
 --no-sex \
 --no-pheno \
 --homozyg-window-snp 50 \
 --homozyg-snp 50 \
 --homozyg-window-missing 3 \
 --homozyg-kb 1 \
 --homozyg-density 1000 

plink \
 --vcf GVCF_SNPs.vcf.gz \
 --recode \
 --out ./plink/GVCF_SNPs \
 --allow-extra-chr  \
 --no-parents \
 --no-sex \
 --no-pheno  \
 --homozyg-window-snp 50 \
 --homozyg-snp 50 \
 --homozyg-window-missing 3 \
 --homozyg-kb 1 \
 --homozyg-density 1000 

plink \
 --vcf GVCF_SNPs.vcf.gz \
 --hardy \
 --out ./plink/GVCF_SNPs \
 --allow-extra-chr  \
 --no-parents \
 --no-sex \
 --no-pheno  \
 --homozyg-window-snp 50 \
 --homozyg-snp 50 \
 --homozyg-window-missing 3 \
 --homozyg-kb 1 \
 --homozyg-density 1000 

