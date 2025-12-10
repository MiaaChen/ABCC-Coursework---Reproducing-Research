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
INPUT_DIR=/scratch/users/k25127820/Group8_ABCC/Cellobiose/alignment_results
OUTPUT_DIR=/scratch/users/k25127820/Group8_ABCC/Cellobiose/Quantify_Gene_Expression
mkdir -p "$OUTPUT_DIR"

REF_GFF=/scratch/users/k25127820/Group8_ABCC/reference_genome/combined_reference/no_fasta_ref.gff
INSERT_GFF=/scratch/users/k25127820/Group8_ABCC/reference_genome/combined_reference/transgenes.gff

# ----- Gene lengths CSV -----
GENE_LENGTHS="$OUTPUT_DIR/gene_lengths.csv"

# ----- Extract gene lengths from GFF -----
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
SAMPLES=(SRR1166445 SRR1166446 SRR1166447)
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
    done
done

echo "HTSeq-count completed for all samples and genomes."

# ----- Normalize CPM + RPKM separately -----
python3 - <<EOF
import pandas as pd
import os

OUTPUT_DIR = "$OUTPUT_DIR"
GENE_LENGTHS = "$GENE_LENGTHS"

samples = ["SRR1166445", "SRR1166446", "SRR1166447"]
genomes = ["ref", "insert"]

lengths_df = pd.read_csv(GENE_LENGTHS, index_col=0)

for s in samples:
    for g in genomes:
        file_path = os.path.join(OUTPUT_DIR, f"{s}_{g}_counts.txt")
        df = pd.read_csv(file_path, sep="\t", header=None)
        df.columns = ["gene", "count"]
        df = df.set_index("gene")

        total_reads = df["count"].sum()
        cpm = df["count"] / total_reads * 1e6

        matched = df.index.intersection(lengths_df.index)
        rpkm = (df.loc[matched, "count"] * 1e9) / (total_reads * lengths_df.loc[matched, "length_bp"])

        df_out = pd.DataFrame({"count": df["count"], "CPM": cpm, "RPKM": rpkm})
        df_out.to_csv(os.path.join(OUTPUT_DIR, f"{s}_{g}_normalized.csv"))

print("CPM and RPKM normalization complete for reference and insert genomes separately.")
EOF

echo "Pipeline finished successfully!"
