#!/bin/bash
#SBATCH --job-name="heterozygosity"
#SBATCH --export=ALL
#SBATCH --mem=300G 


source activate angsd

mkdir heterozygosity
mkdir output

# Put the reference genome in a variable:
REF=$1
#REF="/home/tmichel/projects/rbge/Begonia_genomes/H_sandwicensis_hybrid_assembly_MD_05-2023.fasta"


# SOURCE:
# http://www.popgen.dk/angsd/index.php/Heterozygosity

# List of bam files:
basename -s .sorted.bam -a mapped/*.sorted.bam | sort -V  > mapped/list_bam.txt


# saf files:
while read f ;
do
        angsd \
                -i mapped/"$f".sorted.bam \
                -anc "$REF" \
                -dosaf 1 \
                -gl 1 \
                -out heterozygosity/"$f" ;
        done < mapped/list_bam.txt


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
