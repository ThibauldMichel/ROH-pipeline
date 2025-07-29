#!/bin/bash

#SBATCH --job-name="angsd"
#SBATCH --export=ALL
#SBATCH --mem=200G
#SBATCH --partition=medium


# Load conda environment setup
source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh

# Ensure conda is in PATH for subprocesses (important for snakemake's internal calls)
export PATH=/mnt/apps/users/tmichel/conda/bin:$PATH

# Activate your conda environment
conda activate snakemake_env

# =============================
# CREATE DIRECTORIES
# =============================

mkdir -p heterozygosity
mkdir -p output

# ========================
# CONFIGURATION
# ========================
REF="/home/tmichel/scratch/ROH-pipeline/ROH-pipeline_Bfluvialis_wrapper/Bflu_genome.fa"
COUNTS_DIR="heterozygosity/counts"
SITES_DIR="heterozygosity/sites" 
INTERSECT_DIR="heterozygosity/intersect"
INTERSECT_TMP="${INTERSECT_DIR}/tmp_intersect.bed"
INTERSECT_POSITIONS="${INTERSECT_DIR}/shared_positions.txt"
INTERSECT_BED="${INTERSECT_DIR}/shared_regions.bed"

# Create directories
mkdir -p ${COUNTS_DIR} ${SITES_DIR} ${INTERSECT_DIR}

# =============================
# STEP 0: Create sample list
# =============================
echo "Creating sample list..."
basename -s .sorted.bam -a mapped/*.sorted.bam | sort -V > mapped/sample_list
cp mapped/sample_list mapped/list_bam.txt

# =============================
# STEP 1: Generate callable sites per sample
# =============================
echo "Generating per-sample callable site lists..."

#while read sample; do
#    echo "▶ Processing $sample ..."
#    
#    # Run ANGSD to get callable sites
#    angsd -i mapped/${sample}.sorted.bam \
#        -ref "$REF" \
#        -doCounts 1 \
#        -dumpCounts 2 \
#        -minMapQ 30 \
#        -minQ 20 \
#        -remove_bads 1 \
#        -uniqueOnly 1 \
#        -only_proper_pairs 1 \
#        -out ${COUNTS_DIR}/${sample}
#    
#    # Verify ANGSD ran successfully
#    if [ $? -ne 0 ]; then
#        echo "Error: ANGSD failed for $sample"
#        exit 1
#    fi
#    
#    # Extract positions from .pos.gz file (skip header, convert to 0-based BED)
#    echo "Extracting callable sites..."
#    zcat ${COUNTS_DIR}/${sample}.pos.gz | \
#        awk 'NR>1 {print $1"\t"($2-1)}' | \
#        sort -k1,1 -k2,2n > ${SITES_DIR}/${sample}.sites
#    
#    # Create temporary BED version for intersection
#    awk '{print $1"\t"$2"\t"($2+1)}' ${SITES_DIR}/${sample}.sites > ${SITES_DIR}/${sample}.bed
#    
#done < mapped/sample_list


files=$(cat heterozygosity/intersect/bed_list.txt)


# Start with the first two files
bedtools intersect -a $(echo "$files" | head -n1) -b $(echo "$files" | head -n2 | tail -n1) > temp1.bed

# Sequentially intersect with remaining files
for file in $(echo "$files" | tail -n +3); do
    bedtools intersect -a temp1.bed -b "$file" > temp2.bed
    mv temp2.bed temp1.bed
done

# Final result
mv temp1.bed common_regions.bed

wc -l common_regions.bed
