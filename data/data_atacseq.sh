#!/bin/bash
# download_references.sh
set -e  # Exit on error

echo "Creating directory structure..."
mkdir -p data/reference data/motifs data/02_alignment/bowtie2/index

# Function to check if download was successful
check_download() {
    if [ $? -eq 0 ]; then
        echo "✅ Success: $1"
    else
        echo "❌ Failed: $1"
        exit 1
    fi
}

echo "=== Downloading Reference Genome ==="
# Download human genome (GRCh38) from Ensembl
wget -c https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -O data/reference/genome.fa.gz
check_download "Genome download"
gunzip data/reference/genome.fa.gz
check_download "Genome decompression"

echo "=== Downloading Genome Index for Bowtie2 ==="
# Note: Better to build index locally, but if you want pre-built:
echo "Building Bowtie2 index locally..."
bowtie2-build data/reference/genome.fa data/02_alignment/bowtie2/index/genome
check_download "Bowtie2 index building"

echo "=== Downloading Chromosome Sizes ==="
# Generate chrom.sizes from genome
samtools faidx data/reference/genome.fa
cut -f1,2 data/reference/genome.fa.fai > data/reference/genome.chrom.sizes
check_download "Chromosome sizes generation"

echo "=== Downloading TSS Annotations ==="
# Download Gencode TSS annotations
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz -O data/reference/annotation.gtf.gz
check_download "TSS annotations download"
gunzip data/reference/annotation.gtf.gz
check_download "TSS annotations decompression"

# Create TSS bed file from GTF
awk 'BEGIN{OFS="\t"} $3=="gene" {split($10, g, "\""); print $1, $4-1, $5, g[2], ".", $7}' data/reference/annotation.gtf > data/reference/tss.bed
check_download "TSS bed file creation"

echo "=== Downloading ENCODE Blacklist ==="
wget -c https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz -O data/reference/ENCODE_blacklist.bed.gz
check_download "Blacklist download"
gunzip data/reference/ENCODE_blacklist.bed.gz
check_download "Blacklist decompression"

echo "=== Downloading JASPAR Motifs ==="
mkdir -p data/motifs
wget -c http://jaspar.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt -O data/motifs/jaspar_vertebrates.meme
check_download "JASPAR motifs download"

echo "=== Downloading PhiX Control ==="
wget -c https://downloads.pacbcloud.com/public/references/PhiX174.fa -O data/reference/PhiX174.fa
check_download "PhiX download"

echo "=== Creating Additional Required Files ==="
# Create GFF from GTF (if needed)
cp data/reference/annotation.gtf data/reference/annotation.gff
check_download "GFF file creation"

echo "=== Verification ==="
echo "Checking all files were downloaded:"
ls -la data/reference/
ls -la data/motifs/
ls -la data/02_alignment/bowtie2/index/

echo "=== Summary ==="
echo "✅ All reference files downloaded successfully!"
echo "📍 Location: data/reference/"
echo "📊 Files downloaded:"
echo "   - Genome: Homo_sapiens.GRCh38.dna.primary_assembly.fa"
echo "   - Bowtie2 index: built locally"
echo "   - Chromosome sizes: genome.chrom.sizes"
echo "   - TSS annotations: annotation.gtf and tss.bed"
echo "   - ENCODE blacklist: ENCODE_blacklist.bed"
echo "   - JASPAR motifs: jaspar_vertebrates.meme"
echo "   - PhiX control: PhiX174.fa"

