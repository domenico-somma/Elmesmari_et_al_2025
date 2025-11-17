#!/bin/bash

# renamer_BD.sh - Converts BD Rhapsody BAM file to FASTQ with souporcell-compatible barcodes
#
# This script extracts reads from a BD Rhapsody BAM file and converts them to
# FASTQ format with modified read headers that souporcell can understand.
#
# Input: BD Rhapsody BAM file -> Output: FASTQ file with converted barcodes

# Define input and output files for clarity
INPUT_BAM="Organoids_Day1_Hea.bam"
BARCODE_WHITELIST="BD_barcodes.csv"        # Original numerical BD barcodes
BARCODE_CONVERSION="BD_10X_barcodes.csv"   # Mapping: BD_num -> 10X_nucleotide
OUTPUT_FASTQ="O1.tdo.fq"

# Pipeline steps:
samtools view -h "$INPUT_BAM" | \
  # 1. Filter BAM: keep only reads with barcodes present in the whitelist
  perl /datastore/tdo/cellranger.filterBAMvalidbarcodes.pl "$BARCODE_WHITELIST" | \
  # 2. Convert barcodes and reformat as FASTQ
  perl -e '
    # Load the barcode conversion mapping into a hash
    open F, "'"$BARCODE_CONVERSION"'" or die "file not found\n"; 
    while (<F>){ 
      chomp; 
      /^(\S+)\t(\S+)/; 
      $h{$1}=$2
    }; 
    
    # Process each SAM line from stdin
    while (<STDIN>){ 
      @ar=split(/\t/);
      
      # Extract cell barcode (CB tag) and molecular barcode (MR tag)
      /CB:Z:(\S+)/; 
      $cell=$h{$1}; 
      /MR:Z:(\S+)/; 
      $umi=$1; 
      
      # If both barcodes are valid, output as FASTQ
      if (defined($cell) && defined($umi)){ 
        print "\@$ar[0];$cell;$umi\n$ar[9]\n+\n$ar[10]\n";
      }
    }
  ' > "$OUTPUT_FASTQ"

echo "Conversion complete: $INPUT_BAM -> $OUTPUT_FASTQ"