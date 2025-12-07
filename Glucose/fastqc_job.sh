#!/bin/bash
module load fastqc
mkdir -p fastqc_results

for file in *.fastq *.fastq.gz; do
    if [[ -f "$file" ]]; then
        echo "Running FastQC on $file"
        fastqc "$file" -o fastqc_results -t 4
    fi
done 

