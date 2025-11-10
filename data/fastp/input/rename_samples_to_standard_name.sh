#!/bin/bash

# Navigate to input directory
cd ~/ATACseq_Pipeline/ATACseq_Pipeline/data/fastp/input

# Create mapping: SRR ID → sample name
declare -A mapping=(
    [SRR891270]="sample1"
    [SRR891272]="sample2"
    [SRR891275]="sample3"
    [SRR891277]="sample4"
    [SRR891279]="sample5"
)

# Rename files
for srr in "${!mapping[@]}"; do
    sample="${mapping[$srr]}"
    
    # Rename _1.fastq to _R1.fastq.gz
    if [ -f "${srr}_1.fastq" ]; then
        mv "${srr}_1.fastq" "${sample}_R1.fastq"
        gzip "${sample}_R1.fastq"
        echo "Renamed: ${srr}_1.fastq → ${sample}_R1.fastq.gz"
    fi
    
    # Rename _2.fastq to _R2.fastq.gz
    if [ -f "${srr}_2.fastq" ]; then
        mv "${srr}_2.fastq" "${sample}_R2.fastq"
        gzip "${sample}_R2.fastq"
        echo "Renamed: ${srr}_2.fastq → ${sample}_R2.fastq.gz"
    fi
done

echo "✅ All files renamed!"
ls -lh *.fastq.gz
