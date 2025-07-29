#!/bin/bash

#SBATCH --job-name="jellyfish"
#SBATCH --export=ALL
#SBATCH --mem=200G
#SBATCH --partition=medium


# Load conda environment setup
source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh

# Ensure conda is in PATH for subprocesses (important for snakemake's internal calls)
export PATH=/mnt/apps/users/tmichel/conda/bin:$PATH

# Activate your conda environment
conda activate smudgeplot_env

# Run jellyfish
echo "run jellyfish in a loop with all samples"

mkdir -p jellyfish

#cat mapped/sample_list | while read f;
#do
#	echo "processing jellyfish count "$f"" ;
#	jellyfish count -C  -m 21 -s 10G -t 16 -o jellyfish/kmer_counts_"$f".jf <(zcat reads/"$f"_1.fq.gz) <(zcat reads/"$f"_2.fq.gz) ;
#	echo "processing jellyfish dump "$f"" ;
#      jellyfish dump jellyfish/kmer_counts_"$f".jf > jellyfish/kmers_"$f".fasta ;
#	echo "processing jellyfish histo "$f"" ;
#	jellyfish histo jellyfish/kmer_counts_"$f".jf > jellyfish/kmers_histogram_"$f".txt
#done

# The flags:
#  -C           # Count canonical k-mers (strand-neutral, treats reverse complements as identical)
#  -m 21        # K-mer size = 21 bp (common choice for genome heterozygosity/assembly)
#  -s 10G       # Allocate 10GB of memory for the hash table
#  -t 16        # Use 16 CPU threads
#  -o kmer_counts.jf  # Output file name (binary Jellyfish format)
#  reads_1.fastq reads_2.fastq  # Input paired-end FASTQ files


# Optional: Generate histogram for GenomeScope/Smudgeplot
#jellyfish histo -t 16 kmer_counts.jf > kmer_histogram.txt



echo "make an histogram for genoscope"
#
#cat mapped/sample_list | while read f;
#do
#	jellyfish histo -t 16 jellyfish/kmer_counts_"$f".jf > jellyfish/mer_count_"$f".histo ;
#done
# This gives two columns: k-mer coverage vs. number of k-mers.

echo "create directory genoscope"

#echo "modify samples list for R"
#tr -d '\r' < mapped/sample_list > mapped/sample_list_for_R
#sed -i 's/[ \t]*$//' mapped/sample_list_for_R

echo "compute heterozygosity"
#cat mapped/sample_list_for_R | while read f;
#do
#	mkdir -p jellyfish/genoscope/geno_"$f"
#	Rscript /home/tmichel/genomescope2.0/genomescope.R -i jellyfish/mer_count_"$f".histo -k 21 -o jellyfish/genoscope/geno_"$f" ;
#done

#!/bin/bash

list="mapped/sample_list"
outfile="jellyfish/table_het_jellyfish.csv"

echo "sample,homozygosity,heterozygosity" > "$outfile"

while read sample; do
    # Skip empty lines
    [ -z "$sample" ] && continue

    summary="jellyfish/genoscope/geno_${sample}/summary.txt"

    if [[ -f "$summary" ]]; then
        # Extract min Homozygous and Heterozygous percentages
        homo=$(grep "Homozygous (aa)" "$summary" | awk '{print $3}' | tr -d '%')
        hetero=$(grep "Heterozygous (ab)" "$summary" | awk '{print $4}' | tr -d '%')

        # Convert to ratio (divide by 100)
        homo_ratio=$(awk "BEGIN {printf \"%.6f\", $homo/100}")
        hetero_ratio=$(awk "BEGIN {printf \"%.6f\", $hetero/100}")

        echo "$sample,$homo_ratio,$hetero_ratio" >> "$outfile"
    else
        echo "Warning: summary file missing for $sample" >&2
    fi
done < "$list"

echo "Table saved to $outfile"



