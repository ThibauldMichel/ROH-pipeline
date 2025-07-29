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

# =============================
# VARIABLES
# =============================

# Put the reference genome in a variable:
REF=$1
SITES_DIR="heterozygosity/sites"
COUNTS_DIR="heterozygosity/counts"
INTERSECT_POSITIONS="${SITES_DIR}/sites_intersect.txt"
INTERSECT_BED="${SITES_DIR}/sites_intersect.bed"

mkdir -p $SITES_DIR $COUNTS_DIR

# =============================
# STEP 0: Create sample list
# =============================
basename -s .sorted.bam -a mapped/*.sorted.bam | sort -V > mapped/sample_list
cp mapped/sample_list mapped/list_bam.txt

# =============================
# STEP 1: Generate callable sites per sample
# =============================
echo "Generating per-sample callable site lists..."

while read sample; do
    echo "Processing $sample ..."
    angsd -i mapped/${sample}.sorted.bam \
        -ref "$REF" \
        -dumpCounts 2 \
        -minMapQ 30 \
        -minQ 20 \
        -remove_bads 1 \
        -uniqueOnly 1 \
        -only_proper_pairs 1 \
        -out ${COUNTS_DIR}/${sample}

    # Convert counts to positions (0-based for BED)
    zcat ${COUNTS_DIR}/${sample}.counts.gz | \
        awk '{print $1"\t"($2-1)}' | sort -k1,1 -k2,2n > ${SITES_DIR}/${sample}.sites
done < mapped/sample_list

# =============================
# STEP 2: Intersect callable sites
# =============================
echo "Creating intersection of all callable sites..."

# Start with the first sample
head -n 1 mapped/sample_list | while read first_sample; do
    cp ${SITES_DIR}/${first_sample}.sites $INTERSECT_POSITIONS
done

# Intersect with remaining samples
tail -n +2 mapped/sample_list | while read sample; do
    sort $INTERSECT_POSITIONS > ${SITES_DIR}/tmp_intersect_sorted.txt
    sort ${SITES_DIR}/${sample}.sites > ${SITES_DIR}/tmp_sample_sorted.txt

    comm -12 ${SITES_DIR}/tmp_intersect_sorted.txt ${SITES_DIR}/tmp_sample_sorted.txt > ${SITES_DIR}/tmp_intersect.txt
    mv ${SITES_DIR}/tmp_intersect.txt $INTERSECT_POSITIONS
done

# =============================
# STEP 2b: Convert positions to BED intervals
# =============================
echo "Converting intersected positions to BED intervals..."

awk '
BEGIN {OFS="\t"}
{
  if (NR==1) {
    chr=$1; start=$2; prev=$2
  } else {
    if ($1==chr && $2==prev+1) {
      prev=$2
    } else {
      print chr, start, prev+1
      chr=$1; start=$2; prev=$2
    }
  }
}
END {
  print chr, start, prev+1
}' $INTERSECT_POSITIONS > $INTERSECT_BED

# =============================
# STEP 2c: Index BED file for ANGSD
# =============================
echo "Indexing BED intervals..."
angsd sites index $INTERSECT_BED

# =============================
# STEP 3: Summarize site counts
# =============================
echo -e "Sample\tCallable_Sites\tShared_Sites\tIgnored_Sites" > heterozygosity/summary_sites.txt
shared_sites=$(wc -l < $INTERSECT_POSITIONS)
while read sample; do
    total=$(wc -l < ${SITES_DIR}/${sample}.sites)
    ignored=$((total - shared_sites))
    echo -e "${sample}\t${total}\t${shared_sites}\t${ignored}"
done < mapped/sample_list >> heterozygosity/summary_sites.txt

# =============================
# STEP 4: Generate SAF files using BED
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
done < mapped/sample_list

echo "Done."


# Step 5: Estimate unfolded SFS per sample
echo "Estimating SFS for each sample..."
while read sample; do
    realSFS heterozygosity/${sample}.saf.idx > heterozygosity/${sample}.ml
done < mapped/sample_list

# Step 6: Build output matrix
echo -e "Sample\tHom_Ancestral\tHeterozygous\tHom_Derived\tHeterozygosity" > heterozygosity/matrix_het_angsd.tsv
while read sample; do
    read homA het homD < heterozygosity/${sample}.ml
    total=$(echo "$homA + $het + $homD" | bc -l)
    het_frac=$(echo "$het / $total" | bc -l)
    echo -e "${sample}\t${homA}\t${het}\t${homD}\t${het_frac}"
done < mapped/sample_list >> heterozygosity/matrix_het_angsd.tsv

echo "✔ All steps completed."
