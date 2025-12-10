#!/bin/bash
#SBATCH --job-name=htseq_rpkm
#SBATCH --output=htseq_rpkm.out
#SBATCH --error=htseq_rpkm.err
#SBATCH --time=04:00:00
#SBATCH --mem=16G

# ----- Load modules -----
module load py-htseq/0.11.2-gcc-13.2.0-python-3.11.6
module load py-pandas/1.5.3-gcc-13.2.0-python-3.11.6

# ----- Directories -----
INPUT_DIR=/scratch/users/k25127820/Group8_ABCC/Glucose/alignment_results
OUTPUT_DIR=/scratch/users/k25127820/Group8_ABCC/Glucose/Quantify_Gene_Expression
mkdir -p "$OUTPUT_DIR"

REF_GFF=/scratch/users/k25127820/Group8_ABCC/reference_genome/combined_reference/no_fasta_ref.gff
INSERT_GFF=/scratch/users/k25127820/Group8_ABCC/reference_genome/combined_reference/transgenes.gff

# ----- Gene lengths CSV -----
GENE_LENGTHS="$OUTPUT_DIR/gene_lengths.csv"

# ----- Extract gene lengths from GFF -----
echo "Extracting gene lengths from GFF files..."
python3 - <<EOF
import pandas as pd

# Bash variables
GENE_LENGTHS = "$GENE_LENGTHS"
REF_GFF = "$REF_GFF"
INSERT_GFF = "$INSERT_GFF"

def gff_to_lengths(gff_file, id_attr):
    lengths = []
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            start = int(parts[3])
            end = int(parts[4])
            info = parts[8]
            gene_id = None
            for field in info.split(";"):
                if field.startswith(f"{id_attr}="):
                    gene_id = field.split("=")[1]
                    break
            if gene_id:
                lengths.append((gene_id, end - start + 1))
    return pd.DataFrame(lengths, columns=["gene_id", "length_bp"])

ref_lengths = gff_to_lengths(REF_GFF, "Parent")
insert_lengths = gff_to_lengths(INSERT_GFF, "gene_id")

all_lengths = pd.concat([ref_lengths, insert_lengths])
all_lengths = all_lengths.groupby("gene_id", as_index=False).sum()
all_lengths.to_csv(GENE_LENGTHS, index=False)
print(f"Gene lengths saved to {GENE_LENGTHS}")
EOF

# ----- HTSeq-count for all samples (reference and insert separately) -----
SAMPLES=(SRR1166442 SRR1166443 SRR1166444)
GENOMES=("ref" "insert")

for SAMPLE in "${SAMPLES[@]}"; do
    for GENOME in "${GENOMES[@]}"; do
        if [ "$GENOME" == "ref" ]; then
            GFF="$REF_GFF"
            IDATTR="Parent"
        else
            GFF="$INSERT_GFF"
            IDATTR="gene_id"
        fi

        echo "Counting $SAMPLE for $GENOME genome..."
        htseq-count -f bam -r pos --idattr="$IDATTR" \
            "$INPUT_DIR/${SAMPLE}_sorted.bam" \
            "$GFF" \
            > "$OUTPUT_DIR/${SAMPLE}_${GENOME}_counts.txt"
        
        # Check if output file has content
        if [ ! -s "$OUTPUT_DIR/${SAMPLE}_${GENOME}_counts.txt" ]; then
            echo "WARNING: Empty output for ${SAMPLE}_${GENOME}"
        fi
    done
done

echo "HTSeq-count completed for all samples and genomes."

# ----- Normalize CPM + RPKM separately -----
echo "Starting normalization..."
python3 - <<EOF
import pandas as pd
import os

OUTPUT_DIR = "$OUTPUT_DIR"
GENE_LENGTHS = "$GENE_LENGTHS"

# FIXED: Use the same sample IDs as in HTSeq-count
samples = ["SRR1166442", "SRR1166443", "SRR1166444"]
genomes = ["ref", "insert"]

# Read gene lengths
lengths_df = pd.read_csv(GENE_LENGTHS)
lengths_df = lengths_df.set_index("gene_id")

for s in samples:
    for g in genomes:
        file_path = os.path.join(OUTPUT_DIR, f"{s}_{g}_counts.txt")
        
        # Check if file exists and has content
        if not os.path.exists(file_path):
            print(f"ERROR: File not found: {file_path}")
            continue
        
        if os.path.getsize(file_path) == 0:
            print(f"ERROR: Empty file: {file_path}")
            continue
        
        # Read counts
        df = pd.read_csv(file_path, sep="\t", header=None)
        df.columns = ["gene", "count"]
        df = df.set_index("gene")
        
        # Remove HTSeq special rows
        special_rows = ["__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique"]
        df = df.drop(special_rows, errors='ignore')
        
        # Calculate CPM
        total_reads = df["count"].sum()
        cpm = df["count"] / total_reads * 1e6
        
        # Calculate RPKM (only for genes with known lengths)
        matched = df.index.intersection(lengths_df.index)
        rpkm = pd.Series(index=df.index, dtype=float)
        if len(matched) > 0:
            rpkm.loc[matched] = (df.loc[matched, "count"] * 1e9) / (total_reads * lengths_df.loc[matched, "length_bp"])
        
        # Create output dataframe
        df_out = pd.DataFrame({
            "count": df["count"], 
            "CPM": cpm, 
            "RPKM": rpkm
        })
        
        output_file = os.path.join(OUTPUT_DIR, f"{s}_{g}_normalized.csv")
        df_out.to_csv(output_file)
        print(f"Saved: {output_file}")

print("CPM and RPKM normalization complete for reference and insert genomes separately.")
EOF

echo "Pipeline finished successfully!"
