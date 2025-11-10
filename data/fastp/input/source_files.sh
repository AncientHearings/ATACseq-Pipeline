#!/bin/bash
set -euo pipefail

sra_list=(
    SRR891270
    SRR891272
    SRR891275
    SRR891277
    SRR891279
)

download_sra() {
    sra=$1
    echo "[$(date '+%H:%M:%S')] Starting $sra"
    
    # Step 1: Prefetch from SRA
    if ! prefetch "$sra" -O sra/; then
        echo "[$(date '+%H:%M:%S')] ❌ Prefetch failed for $sra"
        return 1
    fi
    
    # Step 2: Convert to FASTQ (fixed command)
    if ! fasterq-dump "sra/${sra}/${sra}.sra" \
        --split-files \
        --threads 8 \
        --temp /tmp/fasterq-temp-${sra} \
        --outdir fastq/; then
        echo "[$(date '+%H:%M:%S')] ❌ FASTQ dump failed for $sra"
        return 1
    fi
    
    echo "[$(date '+%H:%M:%S')] ✅ Completed $sra"
}

export -f download_sra

# Create directories
mkdir -p fastq sra logs

# Run downloads with 4 parallel jobs
printf "%s\n" "${sra_list[@]}" | xargs -n 1 -P 4 -I {} bash -c 'download_sra "$@"' _ {}

echo "[$(date '+%H:%M:%S')] All downloads completed!"
