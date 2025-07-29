#!/bin/bash

#SBATCH --job-name="snakewrapper"
#SBATCH --export=ALL
#SBATCH --mem=200G
#SBATCH --partition=medium

# Activate conda
#source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh

# Enable conda via the recommended hook
#eval "$(/mnt/apps/users/tmichel/conda/bin/conda shell.bash hook)"


# Activate the conda environment
#conda activate snakemake_env

#echo "version de conda"
#conda --version

#echo "version de snakemake"
#snakemake --version

#echo "which conda"
#which conda

#conda config --set channel_priority strict


#snakemake --cores all --use-conda --latency-wait 60 --rerun-incomplete --conda-cleanup-pkgs cache

#snakemake --cores all --rerun-incomplete --unlock 



# Load conda environment setup
source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh

# Ensure conda is in PATH for subprocesses (important for snakemake's internal calls)
export PATH=/mnt/apps/users/tmichel/conda/bin:$PATH

# Activate your conda environment
conda activate snakemake_env

# re-install mamba
#conda install -c conda-forge mamba=1.5.6 --update-deps --force-reinstall

echo "Snakemake version:"
which snakemake
snakemake --version

echo "Which conda:"
which conda
conda --version

echo "Which mamba"
which mamba
mamba --version

# Run snakemake with all cores, conda support and cleaning up cache
snakemake --cores all --use-conda --latency-wait 60 --rerun-incomplete --conda-cleanup-pkgs cache

#snakemake --cores all --rerun-incomplete --unlock 
