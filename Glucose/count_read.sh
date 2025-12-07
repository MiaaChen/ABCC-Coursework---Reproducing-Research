#!/bin/bash 

#-----load module-----
module load py-htseq/0.11.2-gcc-13.2.0-python-3.11.6  
module load python/3.11.6 
module load py-pandas/1.5.3-gcc-13.2.0-python-3.11.6
#-----directories-----
INPUT_DIR=/scratch/grp/msc_appbio/Group8_ABCC/Glucose/alignment_results
OUTPUT_DIR=/scratch/grp/msc_appbio/Group8_ABCC/Glucose/Quantify_Gene_Expression
REF_GFF=/scratch/grp/msc_appbio/Group8_ABCC/reference_genome/combined_reference/no_fasta_ref.gff
INSERT_GFF=/scratch/grp/msc_appbio/Group8_ABCC/reference_genome/combined_reference/transgenes_clean.gff
MER_INPUT=/scratch/grp/msc_appbio/Group8_ABCC/Cellobiose/Quantify_Gene_Expression

mkdir -p "$OUTPUT_DIR"

echo "Starting HTSeq-count for all samples"

#-----SRR1166445-----
echo "Counting SRR1166442 for reference genome"
htseq-count -f bam -r pos --idattr=Parent \
    "$INPUT_DIR/SRR1166442_sorted.bam" \
    "$REF_GFF" \
    > "$OUTPUT_DIR/SRR1166442_ref_counts.txt"
echo "SRR1166442 reference genome complete"

echo "Counting SRR1166442 for insert genome"
htseq-count -f bam -r pos --idattr=gene_id \
    "$INPUT_DIR/SRR1166442_sorted.bam" \
    "$INSERT_GFF" \
    > "$OUTPUT_DIR/SRR1166442_insert_counts.txt"
echo "SRR1166442 insert genome complete"

#-----merge the counts-----
echo "Merging counts for SRR1166442"

# combine both files
cat "$MER_INPUT/SRR1166442_ref_counts.txt" \
    "$MER_INPUT/SRR1166442_insert_counts.txt" \
    | sort -k1,1 \
    | awk -F'\t' '
        {
            gene=$1;
            count=$2;
            sum[gene] += count;
        }
        END {
            for (g in sum) {
                print g "\t" sum[g];
            }
        }
    ' \
    | sort -k1,1 \
    > "$MER_INPUT/SRR1166442_merged_counts.txt"

echo "SRR1166442 counts merged"


#-----SRR1166443-----
echo "Counting SRR1166443 for reference genome..."
htseq-count -f bam -r pos --idattr=Parent \
    "$INPUT_DIR/SRR1166443_sorted.bam" \
    "$REF_GFF" \
    > "$OUTPUT_DIR/SRR1166443_ref_counts.txt"
echo "SRR1166443 reference genome complete"

echo "Counting SRR1166443 for insert genome..."
htseq-count -f bam -r pos --idattr=gene_id \
    "$INPUT_DIR/SRR1166443_sorted.bam" \
    "$INSERT_GFF" \
    > "$OUTPUT_DIR/SRR1166443_insert_counts.txt"
echo "SRR1166443 insert genome complete"

# Merge counts
echo "Merging counts for SRR1166443"

# combine both files
cat "$MER_INPUT/SRR1166443_ref_counts.txt" \
    "$MER_INPUT/SRR1166443_insert_counts.txt" \
    | sort -k1,1 \
    | awk -F'\t' '
        {
            gene=$1;
            count=$2;
            sum[gene] += count;
        }
        END {
            for (g in sum) {
                print g "\t" sum[g];
            }
        }
    ' \
    | sort -k1,1 \
    > "$MER_INPUT/SRR1166443_merged_counts.txt"

echo "SRR1166443 counts merged"


#-----SRR1166444-----
echo "Counting SRR1166444 for reference genome..."
htseq-count -f bam -r pos --idattr=Parent \
    "$INPUT_DIR/SRR1166444_sorted.bam" \
    "$REF_GFF" \
    > "$OUTPUT_DIR/SRR1166444_ref_counts.txt"
echo "SRR1166444 reference genome complete"

echo "Counting SRR1166444 for insert genome..."
htseq-count -f bam -r pos --idattr=gene_id \
    "$INPUT_DIR/SRR1166444_sorted.bam" \
    "$INSERT_GFF" \
    > "$OUTPUT_DIR/SRR1166444_insert_counts.txt"
echo "SRR1166444 insert genome complete"

# Merge counts
echo "Merging counts for SRR1166445 (Bash only)"

# combine both files
cat "$MER_INPUT/SRR1166444_ref_counts.txt" \
    "$MER_INPUT/SRR1166444_insert_counts.txt" \
    | sort -k1,1 \
    | awk -F'\t' '
        {
            gene=$1;
            count=$2;
            sum[gene] += count;
        }
        END {
            for (g in sum) {
                print g "\t" sum[g];
            }
        }
    ' \
    | sort -k1,1 \
    > "$MER_INPUT/SRR1166444_merged_counts.txt"

echo "SRR1166444 counts merged"
echo "All merged complete"

echo "All HTSeq-count jobs finished!"
