#!/bin/bash

#SBATCH --job-name="snakemake"
#SBATCH --export=ALL
#SBATCH --partition=himem --mem=300G 

#source activate snakemake
#source activate snakemake

source activate /mnt/shared/scratch/tmichel/apps/miniforge3/envs/snakemake

#conda install -c bioconda snakemake-wrapper-utils



snakemake --version





snakemake --cores all --use-conda --latency-wait 60 

#--unlock 
