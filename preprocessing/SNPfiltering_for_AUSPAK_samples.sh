#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J VCF_filter_pipeline_AUSPAKsamples
#SBATCH -o /out/VCF_filter_pipeline_AUSPAKsamples.%J.out
#SBATCH --time=09:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=20

set -e  # Exit on any error

echo "=== Step 1: Extract AUS/PAK samples ==="
echo "Started at $(date)"

module load bcftools

# Subset VCF to samples of interest for genomic precition project
bcftools view CqAllChro_Cp_Mt.SNPs.exclude_filtered_610_samples.vcf.gz \
  -S auspak_sample_list.txt \
  --threads 20 -O z \
  -o CqAllChro_Cp_Mt.SNPs.exclude_filtered_finalAUS_PAK_samples.vcf.gz

echo "Step 1 completed at $(date)"

# Index for downstream processing
bcftools index --threads=20 CqAllChro_Cp_Mt.SNPs.exclude_filtered_finalAUS_PAK_samples.vcf.gz

echo ""
echo "=== Step 2: Filter for quality SNPs ==="
echo "Started at $(date)"

# Apply QC filters for genomic prediction:
# -r: nuclear chromosomes only (exclude organellar)
# -m2 -M2: biallelic SNPs only
# -v snps: SNPs only (no indels)
# F_MISSING < 0.2: max 20% missing data per variant
# MAF >= 0.01: minor allele frequency filter
# MEAN(DP) 5-30: depth filters (avoid low coverage and repeats/paralogs)
bcftools view --threads 20 -r Cq1A,Cq1B,Cq2A,Cq2B,Cq3A,Cq3B,Cq4A,Cq4B,Cq5A,Cq5B,Cq6A,Cq6B,Cq7A,Cq7B,Cq8A,Cq8B,Cq9A,Cq9B \
  -m2 -M2 -v snps \
  -i 'F_MISSING < 0.2 && MAF[0] >= 0.01 && MEAN(FORMAT/DP)>=5 && MEAN(FORMAT/DP)<=30' \
  CqAllChro_Cp_Mt.SNPs.exclude_filtered_finalAUS_PAK_samples.vcf.gz \
  -Oz -o CqAllChro_Cp_Mt.SNPs.exclude_filtered_finalAUS_PAK_samples.chromosomesonly.snpsonly_minmeanDP5-maxmeanDP30_maxmissing0.8_biallelic_MAF0.01.vcf.gz

bcftools index --threads=20 CqAllChro_Cp_Mt.SNPs.exclude_filtered_finalAUS_PAK_samples.chromosomesonly.snpsonly_minmeanDP5-maxmeanDP30_maxmissing0.8_biallelic_MAF0.01.vcf.gz

echo "Step 2 completed at $(date)"

echo ""
echo "=== Step 3: Strip to GT-only format ==="
echo "Started at $(date)"

# Keep only GT field for genomic prediction (removes DP, GQ, etc.)
# Also removes GATK command lines to reduce header clutter
bcftools annotate --threads=20 \
  -x ^FORMAT/GT \
  --header-lines <(zcat CqAllChro_Cp_Mt.SNPs.exclude_filtered_finalAUS_PAK_samples.chromosomesonly.snpsonly_minmeanDP5-maxmeanDP30_maxmissing0.8_biallelic_MAF0.01.vcf.gz | grep "^##" | grep -v "GATKCommandLine") \
  CqAllChro_Cp_Mt.SNPs.exclude_filtered_finalAUS_PAK_samples.chromosomesonly.snpsonly_minmeanDP5-maxmeanDP30_maxmissing0.8_biallelic_MAF0.01.vcf.gz \
  -O v \
  -o CqAllChro_Cp_Mt.SNPs.exclude_filtered_finalAUS_PAK_samples.chromosomesonly.snpsonly_minmeanDP5-maxmeanDP30_maxmissing0.8_biallelic_MAF0.01_GTonly_noGATK.vcf

echo "Step 3 completed at $(date)"
echo ""
echo "=== Pipeline completed successfully! ==="