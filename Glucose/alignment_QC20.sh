#!/bin/bash

#-----------------------------
# Modules
#-----------------------------
module load bowtie2
module load samtools

#-----------------------------
# Paths
#-----------------------------
INDEX=/scratch/grp/msc_appbio/Group8_ABCC/reference_genome/combined_reference/combined_index
INPUT_DIR=/scratch/grp/msc_appbio/Group8_ABCC/Glucose/trimmed_result_threshold_20
OUTPUT_DIR=/scratch/grp/msc_appbio/Group8_ABCC/Glucose/alignment_result_QC20

#-----------------------------
# Threads (default to 4 if not set by SLURM)
#-----------------------------
THREADS=${SLURM_CPUS_PER_TASK:-4}

#-----------------------------
# Create output directory
#-----------------------------
mkdir -p "$OUTPUT_DIR"

#-----------------------------
# Change to input directory
#-----------------------------
cd "$INPUT_DIR" || { echo "ERROR: Cannot change to input directory $INPUT_DIR"; exit 1; }

#-----------------------------
# Loop through FASTQ files
#-----------------------------
for file in *_trimmed.fastq; do
    BASENAME=$(basename "$file" _trimmed.fastq)
    
    echo "-----------------------------------"
    echo "Aligning $file ..."
    
    SAM="$OUTPUT_DIR/${BASENAME}.sam"
    BAM="$OUTPUT_DIR/${BASENAME}.bam"
    SORTED="$OUTPUT_DIR/${BASENAME}_sorted.bam"

    # Bowtie2 alignment
    if bowtie2 -p $THREADS -x $INDEX -U "$file" -S "$SAM"; then
        echo "Bowtie2 alignment finished for $file"

        # Convert SAM -> BAM
        samtools view -@ $THREADS -bS "$SAM" > "$BAM"

        # Sort BAM
        samtools sort -@ $THREADS "$BAM" -o "$SORTED"

        # Index BAM
        samtools index "$SORTED"

        # Remove SAM to save space
        rm "$SAM"

        echo "Done → $SORTED"
    else
        echo "ERROR: Bowtie2 failed for $file. Skipping this sample."
    fi
done

echo "All done!"
