#!/bin/bash

#SBATCH --job-name="repair"
#SBATCH --export=ALL
#SBATCH --mem=200G
#SBATCH --partition=medium

# Activate conda
#source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh

# Enable conda via the recommended hook
#eval "$(/mnt/apps/users/tmichel/conda/bin/conda shell.bash hook)"


# Activate the conda environment
#conda activate snakemake_env


samples=(Bflu_A47 B_pleb Hillebrandia)

for sample in "${samples[@]}"; do
    echo "Repairing sample: $sample"
    repair.sh \
        in1=trimmed/${sample}_1_paired.fastq.gz \
        in2=trimmed/${sample}_2_paired.fastq.gz \
        out1=trimmed/${sample}_1_paired.repair.fastq.gz \
        out2=trimmed/${sample}_2_paired.repair.fastq.gz \
        outs=trimmed/${sample}_singletons.repair.fastq.gz \
        overwrite=true
done
