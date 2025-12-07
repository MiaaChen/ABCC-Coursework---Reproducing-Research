#!/bin/bash 

#-----load module-----
module load py-htseq/0.11.2-gcc-13.2.0-python-3.11.6  
module load python/3.11.6 

#-----directories-----
INPUT_DIR=/scratch/grp/msc_appbio/Group8_ABCC/Cellobiose/alignment_results
OUTPUT_DIR=/scratch/grp/msc_appbio/Group8_ABCC/Cellobiose/Quantify_Gene_Expression
REF_GFF=/scratch/grp/msc_appbio/Group8_ABCC/reference_genome/combined_reference/no_fasta_ref.gff
INSERT_GFF=/scratch/grp/msc_appbio/Group8_ABCC/reference_genome/combined_reference/transgenes_clean.gff

mkdir -p "$OUTPUT_DIR"

echo "Starting HTSeq-count for all samples"

#-----SRR1166445-----
echo "Counting SRR1166445 for reference genome"
htseq-count -f bam -r pos --idattr=Parent \
    "$INPUT_DIR/SRR1166445_sorted.bam" \
    "$REF_GFF" \
    > "$OUTPUT_DIR/SRR1166445_ref_counts.txt"
echo "SRR1166445 reference genome complete"

echo "Counting SRR1166445 for insert genome"
htseq-count -f bam -r pos --idattr=gene_id \
    "$INPUT_DIR/SRR1166445_sorted.bam" \
    "$INSERT_GFF" \
    > "$OUTPUT_DIR/SRR1166445_insert_counts.txt"
echo "SRR1166445 insert genome complete"

#-----merge the counts-----
python3 - <<EOF
import pandas as pd

ref = pd.read_csv("$OUTPUT_DIR/SRR1166445_ref_counts.txt", sep="\t", header=None)
insert = pd.read_csv("$OUTPUT_DIR/SRR1166445_insert_counts.txt", sep="\t", header=None)

ref.columns = ["gene", "count"]
insert.columns = ["gene", "count"]

merged = pd.concat([ref, insert])
merged = merged.groupby("gene", as_index=False).sum()

merged.to_csv("$OUTPUT_DIR/SRR1166445_merged_counts.txt", sep="\t", index=False)
EOF

echo "SRR1166445 counts merged"

#-----SRR1166446-----
echo "Counting SRR1166446 for reference genome..."
htseq-count -f bam -r pos --idattr=Parent \
    "$INPUT_DIR/SRR1166446_sorted.bam" \
    "$REF_GFF" \
    > "$OUTPUT_DIR/SRR1166446_ref_counts.txt"
echo "SRR1166446 reference genome complete"

echo "Counting SRR1166446 for insert genome..."
htseq-count -f bam -r pos --idattr=gene_id \
    "$INPUT_DIR/SRR1166446_sorted.bam" \
    "$INSERT_GFF" \
    > "$OUTPUT_DIR/SRR1166446_insert_counts.txt"
echo "SRR1166446 insert genome complete"

# Merge counts
python3 - <<EOF
import pandas as pd

ref = pd.read_csv("$OUTPUT_DIR/SRR1166446_ref_counts.txt", sep="\t", header=None)
insert = pd.read_csv("$OUTPUT_DIR/SRR1166446_insert_counts.txt", sep="\t", header=None)

ref.columns = ["gene", "count"]
insert.columns = ["gene", "count"]

merged = pd.concat([ref, insert])
merged = merged.groupby("gene", as_index=False).sum()

merged.to_csv("$OUTPUT_DIR/SRR1166446_merged_counts.txt", sep="\t", index=False)
EOF

echo "SRR1166446 counts merged"

#-----SRR1166447-----
echo "Counting SRR1166447 for reference genome..."
htseq-count -f bam -r pos --idattr=Parent \
    "$INPUT_DIR/SRR1166447_sorted.bam" \
    "$REF_GFF" \
    > "$OUTPUT_DIR/SRR1166447_ref_counts.txt"
echo "SRR1166447 reference genome complete"

echo "Counting SRR1166447 for insert genome..."
htseq-count -f bam -r pos --idattr=gene_id \
    "$INPUT_DIR/SRR1166447_sorted.bam" \
    "$INSERT_GFF" \
    > "$OUTPUT_DIR/SRR1166447_insert_counts.txt"
echo "SRR1166447 insert genome complete"

# Merge counts
python3 - <<EOF
import pandas as pd

ref = pd.read_csv("$OUTPUT_DIR/SRR1166447_ref_counts.txt", sep="\t", header=None)
insert = pd.read_csv("$OUTPUT_DIR/SRR1166447_insert_counts.txt", sep="\t", header=None)

ref.columns = ["gene", "count"]
insert.columns = ["gene", "count"]

merged = pd.concat([ref, insert])
merged = merged.groupby("gene", as_index=False).sum()

merged.to_csv("$OUTPUT_DIR/SRR1166447_merged_counts.txt", sep="\t", index=False)
EOF

echo "SRR1166447 counts merged"


echo "All HTSeq-count jobs finished!"
