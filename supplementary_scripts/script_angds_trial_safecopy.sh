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

# =============================
# STEP 2: Intersect callable sites
# =============================
echo "Creating intersection of all callable sites..."

# Initialize with first sample
head -n 1 mapped/sample_list | while read first_sample; do
    cp ${SITES_DIR}/${first_sample}.bed $INTERSECT_TMP
done

# Intersect with remaining samples using bedtools
tail -n +2 mapped/sample_list | while read sample; do
    echo "∩ Intersecting with $sample..."
    bedtools intersect \
        -a $INTERSECT_TMP \
        -b ${SITES_DIR}/${sample}.bed \
        -sorted > ${SITES_DIR}/tmp_intersect.bed
    
    mv ${SITES_DIR}/tmp_intersect.bed $INTERSECT_TMP
done

# Final intersection file
mv $INTERSECT_TMP $INTERSECT_BED

# =============================
# STEP 2b: Merge adjacent positions into BED intervals
# =============================
echo "Merging adjacent positions into BED intervals..."
awk '
BEGIN {OFS="\t"; chr=""; start=0; end=0}
{
    if ($1 != chr) {  # New chromosome
        if (chr != "") print chr, start, end
        chr=$1; start=$2; end=$3
    } else if ($2 > end) {  # Non-adjacent
        print chr, start, end
        chr=$1; start=$2; end=$3
    } else {  # Update end position
        end = $3
    }
}
END {
    if (chr != "") print chr, start, end
}' $INTERSECT_BED > ${INTERSECT_BED}.merged

mv ${INTERSECT_BED}.merged $INTERSECT_BED

# =============================
# STEP 2c: Index BED file for ANGSD
# =============================
echo "Indexing BED file..."
angsd sites index $INTERSECT_BED

# =============================
# STEP 3: Summarize site counts
# =============================
echo "Generating summary statistics..."
shared_sites=$(awk '{sum += $3-$2} END {print sum}' $INTERSECT_BED)

echo -e "Sample\tCallable_Sites\tShared_Sites\tIgnored_Sites" > heterozygosity/summary_sites.txt
while read sample; do
    total=$(wc -l < ${SITES_DIR}/${sample}.sites)
    ignored=$((total - shared_sites))
    echo -e "${sample}\t${total}\t${shared_sites}\t${ignored}"
done < mapped/sample_list >> heterozygosity/summary_sites.txt

# =============================
# STEP 4: Generate SAF files using BED mask
# =============================
echo "Generating SAF files..."

while read sample; do
    echo "▶ Running ANGSD for ${sample}..."
    angsd -i mapped/${sample}.sorted.bam \
        -ref "$REF" \
        -anc "$REF" \
        -GL 2 \
        -doSaf 1 \
        -minMapQ 30 \
        -minQ 20 \
        -remove_bads 1 \
        -uniqueOnly 1 \
        -only_proper_pairs 1 \
        -baq 1 \
        -C 50 \
        -sites "$INTERSECT_BED" \
        -out heterozygosity/${sample}
    
    if [ $? -ne 0 ]; then
        echo "Error: ANGSD SAF generation failed for $sample"
        exit 1
    fi
done < mapped/sample_list

# =============================
# Cleanup temporary files
# =============================
echo "Cleaning up..."
rm -f ${SITES_DIR}/*.bed ${SITES_DIR}/tmp_*

echo "✅ Pipeline completed successfully!"
