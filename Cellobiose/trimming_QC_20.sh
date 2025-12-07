#!/bin/bash
module load trimmomatic 
module load fastqc

# Set directories
INPUT_DIR=/scratch/grp/msc_appbio/Group8_ABCC/Cellobiose
OUTPUT_DIR=/scratch/grp/msc_appbio/Group8_ABCC/Cellobiose/trimmed_result_threshold_20
FASTQC_POST=/scratch/grp/msc_appbio/Group8_ABCC/Cellobiose/fastqc_post_20

# Create necessary directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$FASTQC_POST"

# List of samples
SAMPLES=("SRR1166445" "SRR1166446" "SRR1166447")

# Loop through samples
for SAMPLE in "${SAMPLES[@]}"; do
    INPUT_FILE="$INPUT_DIR/${SAMPLE}.fastq"
    OUTPUT_FILE="$OUTPUT_DIR/${SAMPLE}_trimmed.fastq"

    echo "=========================================="
    echo "Processing $SAMPLE..."

    # Check if input exists
    if [[ ! -r "$INPUT_FILE" ]]; then
        echo "ERROR: Input file $INPUT_FILE not found or not readable. Skipping $SAMPLE."
        continue
    fi

    # Run Trimmomatic SE — **quality trimming only**
    trimmomatic SE -phred33 \
        "$INPUT_FILE" \
        "$OUTPUT_FILE" \
        SLIDINGWINDOW:4:20 \
        MINLEN:30

    # Check if trimming succeeded
    if [[ -f "$OUTPUT_FILE" ]]; then
        echo "$SAMPLE trimming complete: $OUTPUT_FILE"
    else
        echo "ERROR: Trimming failed for $SAMPLE"
    fi
done

echo "=========================================="
echo "Running FastQC on trimmed reads..."
fastqc "$OUTPUT_DIR"/*_trimmed.fastq -o "$FASTQC_POST" -t 3

echo "All trimming and FastQC complete!"
echo "Check FastQC results in: $FASTQC_POST"
