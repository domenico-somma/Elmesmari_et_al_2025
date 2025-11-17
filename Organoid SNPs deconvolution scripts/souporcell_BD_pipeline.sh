#!/bin/bash

# souporcell_BD_pipeline.sh - Full souporcell pipeline for BD Rhapsody data
#
# This script runs the complete souporcell workflow to assign cells to genetic donors.
# It has been adapted to work with BD Rhapsody data that was converted to 10X-compatible format.

# File definitions for the converted BD Rhapsody data

BD_10X=BD_10X_barcodes.csv    # Mapping: BD numerical barcode -> 10X nucleotide sequence
BD=BD_barcodes.csv            # Original BD numerical barcodes
b10X=10X_barcodes.csv         # 10X-compatible nucleotide barcodes (for souporcell input)
FQ=O1.tdo.fq                  # Input FASTQ file with converted barcodes
Genome=/datastore/Dom/souporcell_datasets/files/GRCh38.primary_assembly.genome.fa
N_clusters=5                  # Number of expected donors (clusters)

# -------------------------------------------------------------------
# STAGE 1: Read Alignment
# -------------------------------------------------------------------
# Map reads to reference genome using minimap2 with splicing-aware parameters
# This aligns the converted FASTQ reads to the human genome
minimap2 -ax splice -t 8 -G50k -k 21 -w 11 --sr -A2 -B8 -O12,32 -E2,1 -r200 -p.5 -N20 -f1000,5000 -n2 -m20 -s40 -g2000 -2K50m --secondary=no $Genome $FQ -o minimap.sam

# -------------------------------------------------------------------
# STAGE 2: BAM File Processing
# -------------------------------------------------------------------
# Retag the SAM file to prepare for souporcell (adds necessary tags)
retag.py --sam minimap.sam --out minitagged.bam

# Sort and index the BAM file for downstream processing
samtools sort minitagged.bam -o minitagged_sorted.bam
samtools index minitagged_sorted.bam

# -------------------------------------------------------------------
# STAGE 3: Variant Calling
# -------------------------------------------------------------------
# Call genetic variants from the aligned reads using freebayes
# This identifies SNPs that can distinguish between different donors
freebayes -f $Genome -iXu -C 2 -q 5 -E 1 -m 1 --min-coverage 6 minitagged_sorted.bam > freebayes.vcf

# -------------------------------------------------------------------
# STAGE 4: Souporcell Analysis
# -------------------------------------------------------------------
# Create allele count matrices for reference and alternative alleles
# vartrix counts how many reads support ref/alt alleles at each variant position
vartrix --umi -b minitagged_sorted.bam -c $b10X --scoring-method coverage --threads 5 --ref-matrix ref.mtx --out-matrix alt.mtx -v freebayes.vcf --fasta $Genome

# Run souporcell to cluster cells by genetic identity
# Uses the allele count matrices to assign cells to N genetic donors
souporcell -a alt.mtx -r ref.mtx -b $b10X -k $N_clusters -t 3 > clusters_tmp.tsv 2>clusters_tmp.tsv.err

# Refine the clustering with troublet (souporcell's troubleshooting tool)
troublet -a alt.mtx -r ref.mtx --clusters clusters_tmp.tsv > clusters_10X_barcode.tsv 2>clusters.tsv.err

# -------------------------------------------------------------------
# STAGE 5: Final Output Processing
# -------------------------------------------------------------------
# Add header to BD barcodes file for merging
{ echo 'barcode'; cat $BD; } > BD_barcodes_moment.csv

# Merge the BD numerical barcodes with their cluster assignments
# Creates final output: BD_barcode -> cluster assignment
paste BD_barcodes_moment.csv  clusters_10X_barcode.tsv | column -s '\t' > clusters.tsv

# Extract only the high-confidence singlets (cells assigned to one donor)
awk -F '\t' ' $2 == "singlet" { print $0 }' clusters.tsv > clusters_singlets.tsv
