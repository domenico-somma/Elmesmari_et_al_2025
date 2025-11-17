#!/bin/bash

# SNP Concordance Analysis for Genetic Donor Matching
#
# PURPOSE: Validate souporcell donor assignments by checking genetic concordance
# between clusters from different organoid experiments. This confirms whether
# the same genetic donors appear across multiple experiments.
#
# METHOD: Use Picard's GenotypeConcordance to compare SNP profiles between
# clusters identified in different organoid samples.

# Activate virtual environment with required tools
source ~/anaconda3/bin/activate javaSNP

# -------------------------------------------------------------------
# STAGE 1: Quality Filtering of SNP Calls
# -------------------------------------------------------------------
# Filter VCF files to retain only high-quality SNPs (QUAL > 30)
# This ensures reliable genotype comparisons by removing low-confidence calls
bcftools filter -O v -o O1_qual30.vcf -i '%QUAL>30' O1_cluster_genotypes.vcf
bcftools filter -O v -o O2_qual30.vcf -i '%QUAL>30' O2_cluster_genotypes.vcf

# -------------------------------------------------------------------
# STAGE 2: Cross-Experiment Donor Concordance Analysis
# -------------------------------------------------------------------
# Compare each cluster from Organoid1 with each cluster from Organoid2
# to identify matching genetic donors across experiments
#
# Input files:
# - O1_qual30.vcf: Filtered SNPs from Organoid1 experiment
# - O2_qual30.vcf: Filtered SNPs from Organoid2 experiment
#
# Parameters:
# - TS = Test Sample index (sample in the first VCF file)
# - CS = Comparison Sample index (sample in the second VCF file)
# - O = Output directory/file prefix
#
# This generates a matrix of all possible cluster comparisons to find
# genetic matches between experiments

echo "Running genotype concordance analysis..."

# Organoid1 vs Organoid2 sample comparison
# Compare each cluster from O1 (TS=0-4) with each cluster from O2 (CS=0-3)
# This identifies which genetic donors are shared between experiments

# O1 Cluster 0 vs all O2 clusters
picard GenotypeConcordance TV=O1_qual30.vcf CV=O2_qual30.vcf O=O1vsO2/0vs0 TS=0 CS=0
picard GenotypeConcordance TV=O1_qual30.vcf CV=O2_qual30.vcf O=O1vsO2/0vs1 TS=0 CS=1
picard GenotypeConcordance TV=O1_qual30.vcf CV=O2_qual30.vcf O=O1vsO2/0vs2 TS=0 CS=2
picard GenotypeConcordance TV=O1_qual30.vcf CV=O2_qual30.vcf O=O1vsO2/0vs3 TS=0 CS=3

# O1 Cluster 1 vs all O2 clusters
picard GenotypeConcordance TV=O1_qual30.vcf CV=O2_qual30.vcf O=O1vsO2/1vs0 TS=1 CS=0
picard GenotypeConcordance TV=O1_qual30.vcf CV=O2_qual30.vcf O=O1vsO2/1vs1 TS=1 CS=1
picard GenotypeConcordance TV=O1_qual30.vcf CV=O2_qual30.vcf O=O1vsO2/1vs2 TS=1 CS=2
picard GenotypeConcordance TV=O1_qual30.vcf CV=O2_qual30.vcf O=O1vsO2/1vs3 TS=1 CS=3
picard GenotypeConcordance TV=O1_qual30.vcf CV=O2_qual30.vcf O=O1vsO2/1vs4 TS=1 CS=4

# O1 Cluster 2 vs all O2 clusters
picard GenotypeConcordance TV=O1_qual30.vcf CV=O2_qual30.vcf O=O1vsO2/2vs0 TS=2 CS=0
picard GenotypeConcordance TV=O1_qual30.vcf CV=O2_qual30.vcf O=O1vsO2/2vs1 TS=2 CS=1
picard GenotypeConcordance TV=O1_qual30.vcf CV=O2_qual30.vcf O=O1vsO2/2vs2 TS=2 CS=2
picard GenotypeConcordance TV=O1_qual30.vcf CV=O2_qual30.vcf O=O1vsO2/2vs3 TS=2 CS=3

# O1 Cluster 3 vs all O2 clusters
picard GenotypeConcordance TV=O1_qual30.vcf CV=O2_qual30.vcf O=O1vsO2/3vs0 TS=3 CS=0
picard GenotypeConcordance TV=O1_qual30.vcf CV=O2_qual30.vcf O=O1vsO2/3vs1 TS=3 CS=1
picard GenotypeConcordance TV=O1_qual30.vcf CV=O2_qual30.vcf O=O1vsO2/3vs2 TS=3 CS=2
picard GenotypeConcordance TV=O1_qual30.vcf CV=O2_qual30.vcf O=O1vsO2/3vs3 TS=3 CS=3

# O1 Cluster 4 vs all O2 clusters
picard GenotypeConcordance TV=O1_qual30.vcf CV=O2_qual30.vcf O=O1vsO2/4vs0 TS=4 CS=0
picard GenotypeConcordance TV=O1_qual30.vcf CV=O2_qual30.vcf O=O1vsO2/4vs1 TS=4 CS=1
picard GenotypeConcordance TV=O1_qual30.vcf CV=O2_qual30.vcf O=O1vsO2/4vs2 TS=4 CS=2
picard GenotypeConcordance TV=O1_qual30.vcf CV=O2_qual30.vcf O=O1vsO2/4vs3 TS=4 CS=3

echo "Concordance analysis complete. Results in O1vsO2/ directory."
echo "High concordance scores indicate the same genetic donor across experiments."