#!/bin/bash

# Download example reference files for testing

set -e

echo "Setting up test data and reference files..."

# Create directories
mkdir -p resources/genome
mkdir -p resources/adapters
mkdir -p data

# Download adapter sequences
echo "Downloading adapter sequences..."
wget -O resources/adapters/TruSeq3-PE-2.fa \
    https://github.com/usadellab/Trimmomatic/raw/main/adapters/TruSeq3-PE-2.fa

# Download mouse mm10 reference (chromosome 19 only for testing)
echo "Downloading mouse mm10 chromosome 19 for testing..."
cd resources/genome

# Download chr19 fasta
wget -O chr19.fa.gz \
    http://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chr19.fa.gz
gunzip chr19.fa.gz
mv chr19.fa genome.fa

# Download GTF for chr19
wget -O gencode.vM25.annotation.gtf.gz \
    http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf.gz
gunzip gencode.vM25.annotation.gtf.gz

# Filter GTF for chr19 only
grep "^chr19" gencode.vM25.chr_patch_hapl_scaff.annotation.gtf > genome.gtf
rm gencode.vM25.chr_patch_hapl_scaff.annotation.gtf

# Create TSS bed file
echo "Creating TSS bed file..."
awk '$3=="gene" {
    if($7=="+") {
        print $1"\t"$4-1"\t"$4"\t"$10"\t"$14"\t"$7
    } else {
        print $1"\t"$5-1"\t"$5"\t"$10"\t"$14"\t"$7
    }
}' genome.gtf | sed 's/[";]//g' | sort -k1,1 -k2,2n > tss.bed

cd ../..

echo "Test reference files downloaded successfully!"
echo "Note: This uses only chromosome 19 for testing purposes."
echo "For full analysis, download the complete reference genome."

# Update config for test data
echo "Updating config for test data..."
sed -i.bak 's/effective_size: 2700000000/effective_size: 61431566/' config/config.yaml

echo "Setup complete! You can now run the pipeline with test data."