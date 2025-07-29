#!/bin/bash
#SBATCH --job-name="install numpy"
#SBATCH --export=ALL
#SBATCH --mem=200G
#SBATCH --partition=medium

# ===== 1. SETUP ENVIRONMENT =====
source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh
export PATH=/mnt/apps/users/tmichel/conda/bin:$PATH

# Activate environment and check dependencies
conda activate smudgeplot_env || {
    echo "ERROR: Failed to activate smudgeplot_env"
    exit 1
}

# Check for numpy (required by smudgeplot)
python -c "import numpy" 2>/dev/null || {
    echo "Installing missing dependencies..."
    conda install -y numpy matplotlib scipy || {
        echo "ERROR: Failed to install dependencies"
        exit 1
    }
}



# ===== 2. RUN JELLYFISH =====
mkdir -p jellyfish
OUTPUT_CSV="jellyfish/table_heterozygosity.csv"
echo "sample,heterozygosity,ploidy" > "$OUTPUT_CSV"

while read f; do
    # ===== 3. RUN SMUDGEPLOT =====
    # --- hetkmers ---
    echo "Running smudgeplot hetkmers..."
    smudgeplot.py hetkmers -o "smudgeplot_out_$f" "jellyfish/kmers_$f.fasta" 2> "jellyfish/${f}_hetkmers.log" || {
        echo "ERROR: smudgeplot hetkmers failed for $f"
        echo "$f,ERROR_HETKMERS,NA" >> "$OUTPUT_CSV"
        continue
    }

    # --- plot ---
    echo "Running smudgeplot plot..."
    smudgeplot.py plot -o "smudgeplot_$f" "smudgeplot_out_${f}_coverages.tsv" 2> "jellyfish/${f}_plot.log" || {
        echo "ERROR: smudgeplot plot failed for $f"
        echo "$f,ERROR_PLOT,NA" >> "$OUTPUT_CSV"
        continue
    }

    # ===== 4. EXTRACT RESULTS =====
    SUMMARY_FILE="smudgeplot_${f}_summary.tsv"
    if [[ -f "$SUMMARY_FILE" ]]; then
        heterozygosity=$(awk -F'\t' 'NR==2 {print $2}' "$SUMMARY_FILE")
        ploidy=$(awk -F'\t' 'NR==2 {print $3}' "$SUMMARY_FILE")
        echo "$f,$heterozygosity,$ploidy" >> "$OUTPUT_CSV"
        echo "SUCCESS: $f - heterozygosity=$heterozygosity, ploidy=$ploidy"
    else
        echo "WARNING: No summary file for $f"
        echo "$f,NO_SUMMARY,NA" >> "$OUTPUT_CSV"
    fi
done < mapped/sample_list

echo "==== Final results ===="
cat "$OUTPUT_CSV"
