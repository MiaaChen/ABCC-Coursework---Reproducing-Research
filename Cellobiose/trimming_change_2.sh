#!/bin/bash

module load trimmomatic
module load fastqc

# Set directories
INPUT_DIR=/scratch/grp/msc_appbio/Group8_ABCC/Cellobiose
OUTPUT_DIR=/scratch/grp/msc_appbio/Group8_ABCC/Cellobiose/trimmed_changed_results_2
FASTQC_POST=/scratch/grp/msc_appbio/Group8_ABCC/Cellobiose/fastqc_post_changed_2

mkdir -p "$OUTPUT_DIR"
mkdir -p "$FASTQC_POST"

echo "Starting trimming for all samples"

###############################################
# SRR1166445
###############################################
echo "Trimming SRR1166445..."
trimmomatic SE -phred33 \
$INPUT_DIR/SRR1166445.fastq \
$OUTPUT_DIR/SRR1166445_trimmed.fastq \
HEADCROP:4 \
CROP:31 \
LEADING:5 \
TRAILING:5 \
SLIDINGWINDOW:4:20 \
MINLEN:25

echo "SRR1166445 complete."


###############################################
# SRR1166446
###############################################
echo "Trimming SRR1166446..."
trimmomatic SE -phred33 \
$INPUT_DIR/SRR1166446.fastq \
$OUTPUT_DIR/SRR1166446_trimmed.fastq \
HEADCROP:4 \
CROP:35 \
LEADING:5 \
TRAILING:5 \
SLIDINGWINDOW:4:20 \
MINLEN:25

echo "SRR1166446 complete."


###############################################
# SRR1166447
###############################################
echo "Trimming SRR1166447..."
trimmomatic SE -phred33 \
$INPUT_DIR/SRR1166447.fastq \
$OUTPUT_DIR/SRR1166447_trimmed.fastq \
HEADCROP:4 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:20 \
MINLEN:30

echo "SRR1166447 complete"


###############################################
# FASTQC
###############################################
echo "All trimming complete! Running FastQC on trimmed reads..."
fastqc $OUTPUT_DIR/*_trimmed.fastq -o $FASTQC_POST -t 3

echo "=========================================="
echo "Trimming Summary:"
echo "SRR1166445: Trimmed first 4bp, kept bp 5-35 (~31bp reads)"
echo "SRR1166446: Trimmed first 4bp, kept bp 5-35 (~31bp reads)"
echo "SRR1166447: Trimmed first 4bp, light quality filtering (~46bp reads)"
echo "=========================================="
echo "Check FastQC results in: $FASTQC_POST"
