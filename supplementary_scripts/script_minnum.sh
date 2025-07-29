#!/bin/bash

#SBATCH --job-name="snakewrapper"
#SBATCH --export=ALL
#SBATCH --mem=200G
#SBATCH --partition=medium


# Load conda environment setup
source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh

# Ensure conda is in PATH for subprocesses (important for snakemake's internal calls)
export PATH=/mnt/apps/users/tmichel/conda/bin:$PATH

# Activate your conda environment
conda activate snakemake_env

echo "Snakemake version:"
which snakemake
snakemake --version

echo "Which conda:"
which conda
conda --version

echo "Which mamba"
which mamba
mamba --version

python3 workflow/script_minimum_number_SNPs_for_ROH.py

