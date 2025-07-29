#!/bin/bash
#SBATCH --job-name="het exp angsd"
#SBATCH --export=ALL
#SBATCH --mem=300G 


source activate angsd

mkdir heterozygosity
mkdir output

# Put the reference genome in a variable:
#REF=$1
REF="/home/tmichel/projects/rbge/Begonia_genomes/H_sandwicensis_hybrid_assembly_MD_05-2023.fasta"


# SOURCE:
# http://www.popgen.dk/angsd/index.php/Heterozygosity

# List of bam files:
ls mapped/*.sorted.bam | sort -V  > mapped/bamlist.txt

https://www.popgen.dk/angsd/index.php/Heterozygosity
# 1. Generate SFS
#angsd -b mapped/bamlist.txt -doSaf 1 -anc "$REF" -GL 1 -out heterozygosity/output.saf


# 2. Estimate the SFS:
realSFS heterozygosity/output.saf.idx -fold 1 > heterozygosity/output.ml

# 3. Calculate hererozygosity
#angsd -b mapped/bamlist.txt -doThetas 1 -doSaf 1 -pest heterozygosity/output.sfs -anc "$REF" -GL 1 -out heterozygosity/output

# 4. Extract heterozygosity estimates
#thetaStat do_stat heterozygosity/output.thetas.idx









# Make list files:
basename -s .sorted.bam mapped/*.sorted.bam  | sort -V  > mapped/sample_list

# heterozygosity content:
while read f ;
do
        realSFS -nSites 1000000   heterozygosity/"$f".saf.idx > heterozygosity/"$f".ml ;
done < mapped/sample_list
        
        
        
# Write the zygosity matrix:    
while read f ;
do      
        echo "$f" | tr "\n" "\t" ;
        cat heterozygosity/"$f".ml 
done > heterozygosity/matrix_het_angsd.tsv  < mapped/sample_list
