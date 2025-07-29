#!/bin/bash

# ==============================================================================
# FastQC Adapter Check Script
# Runs FastQC on all FASTQ files and checks for adapter contamination.
# Outputs a summary report for quick diagnostics.
# ==============================================================================

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------
echo "🔧 Setting up directories..."
mkdir -p fastqc_raw  # For FastQC output
mkdir -p adapter_reports  # For extracted adapter summaries

# ------------------------------------------------------------------------------
# 1. Run FastQC on All Samples
# ------------------------------------------------------------------------------
echo "🚀 Running FastQC on all FASTQ files..."
fastqc -o fastqc_raw  *.fq.gz  # Adjust threads (`-t`) as needed

# ------------------------------------------------------------------------------
# 2. Check for Adapter Contamination
# ------------------------------------------------------------------------------
echo "🔍 Checking for adapter contamination..."
for zip_file in fastqc_raw/*_fastqc.zip; do
    sample=$(basename "$zip_file" | sed 's/_fastqc.zip//')
    unzip -q "$zip_file" -d fastqc_raw/  # Extract FastQC report
    
    # Parse adapter content from FastQC data
    grep -A 10 "Adapter Content" "fastqc_raw/${sample}_fastqc/fastqc_data.txt" \
        > "adapter_reports/${sample}_adapter_report.txt"
    
    # Summarize adapter contamination (any adapter >5% is flagged)
    awk '/Illumina Universal Adapter/{if ($2 > 5) print "  ⚠️  " $1 " adapter: " $2 "%"}' \
        "adapter_reports/${sample}_adapter_report.txt" \
        >> "adapter_reports/summary.txt"
done

# ------------------------------------------------------------------------------
# 3. Generate Final Report
# ------------------------------------------------------------------------------
echo "📊 Generating summary report..."
echo "======================================================"
echo "📑 FastQC Adapter Contamination Summary"
echo "======================================================"
echo "Samples with adapter contamination >5%:"
cat "adapter_reports/summary.txt" 2>/dev/null || echo "  ✅ No significant adapter contamination detected."
echo "------------------------------------------------------"
echo "🔗 Detailed reports:"
echo "- FastQC HTML files: fastqc_raw/*_fastqc.html"
echo "- Adapter reports: adapter_reports/*_adapter_report.txt"
echo "======================================================"

# Clean up extracted FastQC folders (optional)
rm -rf fastqc_raw/*_fastqc
