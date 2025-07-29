#!/bin/bash
#SBATCH --job-name=check_paired_reads
#SBATCH --output=logs/check_paired_reads_%j.log
#SBATCH --error=logs/check_paired_reads_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

# Load any modules if needed, e.g.
# module load gzip

samples=(Bflu_A47 B_pleb Hillebrandia)

echo "Starting read pairing check at $(date)"
for sample in "${samples[@]}"; do
    echo "Sample: $sample"

    # Extract read names (header lines) from R1 and R2
    zcat "trimmed/${sample}_1_paired.fastq.gz" | sed -n '1~4p' > /tmp/${sample}_1_names.txt
    zcat "trimmed/${sample}_2_paired.fastq.gz" | sed -n '1~4p' > /tmp/${sample}_2_names.txt

    # Sort read names for comparison
    sort /tmp/${sample}_1_names.txt > /tmp/${sample}_1_sorted.txt
    sort /tmp/${sample}_2_names.txt > /tmp/${sample}_2_sorted.txt

    # Count total reads in each file
    total_r1=$(wc -l < /tmp/${sample}_1_sorted.txt)
    total_r2=$(wc -l < /tmp/${sample}_2_sorted.txt)

    # Count reads common in both files (paired)
    paired=$(comm -12 /tmp/${sample}_1_sorted.txt /tmp/${sample}_2_sorted.txt | wc -l)

    # Reads unique to R1 (unpaired)
    unpaired_r1=$(comm -23 /tmp/${sample}_1_sorted.txt /tmp/${sample}_2_sorted.txt | wc -l)

    # Reads unique to R2 (unpaired)
    unpaired_r2=$(comm -13 /tmp/${sample}_1_sorted.txt /tmp/${sample}_2_sorted.txt | wc -l)

    echo "  Total reads in R1: $total_r1"
    echo "  Total reads in R2: $total_r2"
    echo "  Paired reads: $paired"
    echo "  Unpaired reads in R1: $unpaired_r1"
    echo "  Unpaired reads in R2: $unpaired_r2"
    echo ""

    # Cleanup temporary files
    rm /tmp/${sample}_1_names.txt /tmp/${sample}_2_names.txt /tmp/${sample}_1_sorted.txt /tmp/${sample}_2_sorted.txt
done
echo "Finished read pairing check at $(date)"

