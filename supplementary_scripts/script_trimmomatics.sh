#!/bin/bash

# Input files
R1=./reads/Bflu_A39_1.fq.gz
R2=./reads/Bflu_A39_2.fq.gz

# Output directory
OUTDIR=./trimmed
mkdir -p "$OUTDIR"

# Adapter file bundled with Trimmomatic
ADAPTERS=/home/tmichel/Trimmomatic-0.39/adapters/TruSeq3-PE.fa

# Output files
P1=$OUTDIR/Bflu_A39_1_paired.fastq.gz
U1=$OUTDIR/Bflu_A39_1_unpaired.fastq.gz
P2=$OUTDIR/Bflu_A39_2_paired.fastq.gz
U2=$OUTDIR/Bflu_A39_2_unpaired.fastq.gz

# Run Trimmomatic
trimmomatic PE -threads 4 -phred33 \
  "$R1" "$R2" \
  "$P1" "$U1" "$P2" "$U2" \
  ILLUMINACLIP:"$ADAPTERS":2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

