#!/bin/bash

module load trimmomatic
module load fastqc

# Set directories
INPUT_DIR=/scratch/grp/msc_appbio/Group8_ABCC/Glucose
OUTPUT_DIR=/scratch/grp/msc_appbio/Group8_ABCC/Glucose/trimmed_changed_results_2
FASTQC_POST=/scratch/grp/msc_appbio/Group8_ABCC/Glucose/fastqc_post_changed_2

mkdir -p "$OUTPUT_DIR"
mkdir -p "$FASTQC_POST"

echo "Starting trimming for all samples"

###############################################
# SRR1166445
###############################################
echo "Trimming SRR1166442..."
trimmomatic SE -phred33 \
$INPUT_DIR/SRR1166442.fastq \
$OUTPUT_DIR/SRR1166442_trimmed.fastq \
HEADCROP:4 \
CROP:31 \
LEADING:5 \
TRAILING:5 \
SLIDINGWINDOW:4:20 \
MINLEN:25

echo "SRR1166442 complete."


###############################################
# SRR1166446
###############################################
echo "Trimming SRR1166443..."
trimmomatic SE -phred33 \
$INPUT_DIR/SRR1166443.fastq \
$OUTPUT_DIR/SRR1166443_trimmed.fastq \
HEADCROP:4 \
CROP:38 \
LEADING:5 \
TRAILING:5 \
SLIDINGWINDOW:4:20 \
MINLEN:25

echo "SRR1166443 complete."


###############################################
# SRR1166447
###############################################
echo "Trimming SRR1166444..."
trimmomatic SE -phred33 \
$INPUT_DIR/SRR1166444.fastq \
$OUTPUT_DIR/SRR1166444_trimmed.fastq \
HEADCROP:4 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:20 \
MINLEN:30

echo "SRR1166444 complete"


###############################################
# FASTQC
###############################################
echo "All trimming complete! Running FastQC on trimmed reads..."
fastqc $OUTPUT_DIR/*_trimmed.fastq -o $FASTQC_POST -t 3

echo "=========================================="
echo "Trimming Summary:"
echo "SRR1166442: Trimmed first 4bp, kept bp 5-35 (~31bp reads)"
echo "SRR1166443: Trimmed first 4bp, kept bp 5-35 (~31bp reads)"
echo "SRR1166444: Trimmed first 4bp, light quality filtering (~46bp reads)"
echo "=========================================="
echo "Check FastQC results in: $FASTQC_POST"
