#!/bin/bash

# ------------------------------------------------------------------------------
# Comprehensive FASTQ Quality Check Script
# Checks: Paired-end integrity, read counts, read lengths, and trimming clues.
# Usage: ./check_fastq_metrics.sh
# ------------------------------------------------------------------------------

echo "======================================================"
echo "🔍 FASTQ File Analysis Report"
echo "======================================================"
echo "📌 Note:"
echo "- Paired-end files should have matching _1/_2 reads."
echo "- Trimmed reads often have uniform lengths and no adapters."
echo "- Adapter contamination suggests untrimmed data."
echo "------------------------------------------------------"

# ==============================================================================
# 1. Check Paired-End Integrity
# ==============================================================================
echo "🔄 Checking paired-end file consistency..."
paired_errors=0
for file in *_1.fq.gz; do
    sample=${file%_1.fq.gz}
    if [[ ! -f "${sample}_2.fq.gz" ]]; then
        echo "  ❌ ERROR: ${sample}_2.fq.gz missing (single-end?)"
        paired_errors=$((paired_errors + 1))
    fi
done
if [[ $paired_errors -eq 0 ]]; then
    echo "  ✅ All samples have paired _1/_2 files."
else
    echo "  ⚠️  WARNING: $paired_errors samples lack paired files!"
fi

# ==============================================================================
# 2. Count Reads per File
# ==============================================================================
echo "📊 Counting reads per file (4 lines = 1 read)..."
declare -A read_counts
for file in *.fq.gz; do
    count=$(zgrep -c "^@" "$file")
    read_counts[$file]=$count
    echo "  ${file}: $count reads"
done

# ==============================================================================
# 3. Verify Read-Pair Matching
# ==============================================================================
echo "🔗 Checking if _1/_2 read counts match..."
mismatches=0
for sample in $(ls *_1.fq.gz | sed 's/_1.fq.gz//'); do
    count1=${read_counts["${sample}_1.fq.gz"]}
    count2=${read_counts["${sample}_2.fq.gz"]}
    if [[ $count1 -ne $count2 ]]; then
        echo "  ❌ ${sample}: UNBALANCED (${count1} vs ${count2} reads)"
        mismatches=$((mismatches + 1))
    else
        echo "  ✅ ${sample}: Balanced (${count1} reads each)"
    fi
done
if [[ $mismatches -gt 0 ]]; then
    echo "  ⚠️  WARNING: $mismatches samples have mismatched pairs!"
fi

# ==============================================================================
# 4. Check Read Lengths (Trimming Clues)
# ==============================================================================
echo "📏 Checking read lengths (uniform lengths suggest trimming)..."
for file in *_1.fq.gz; do
    sample=${file%_1.fq.gz}
    echo "  📌 Sample: $sample"
    echo "    Read lengths (top 5):"
    zcat "$file" | awk 'NR%4==2 {print length($0)}' | sort | uniq -c | sort -nr | head -n 5
done

# ==============================================================================
# 5. Adapter Check (Quick Scan)
# ==============================================================================
echo "🧬 Scanning for adapter sequences (first 50 reads)..."
common_adapters=("AGATCGGAAGAGC" "CTGTCTCTTATAC" "GTGACTGGAGTTC")
for file in *_1.fq.gz; do
    echo "  🔍 ${file}:"
    zcat "$file" | head -n 200 | grep -E "${common_adapters[*]}" | \
        awk '{print "    Potential adapter: " $0}'
done | grep -v "Potential adapter" || echo "    No obvious adapters detected."

# ==============================================================================
# Final Report
# ==============================================================================
echo "------------------------------------------------------"
echo "📋 Summary:"
echo "- Paired-end integrity: $paired_errors errors"
echo "- Read-pair matching: $mismatches mismatches"
echo "- Read lengths: Check above for uniformity (trimming clue)."
echo "- Adapters: Manually verify if any were detected."
echo "======================================================"
