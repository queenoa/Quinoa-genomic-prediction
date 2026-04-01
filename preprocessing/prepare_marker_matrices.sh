#!/bin/bash
# prepare_marker_matrices.sh
#
# Export marker data from the quality-filtered VCF into PLINK .raw format
# for the genomic prediction models:
#   1. Full marker set for RKHS (all 1,824,377 SNPs)
#   2. LD-pruned marker set for BayesC (~557k SNPs, r2 < 0.5)
#   3. Small test subset (1000 markers) for unit/integration tests
#
# Input:  quinoa_551accessions_genomic_prediction.vcf (from SNPfiltering)
# Requires: plink2

VCF="quinoa_551accessions_genomic_prediction.vcf"

# Assign unique variant IDs (chr:pos:ref:alt) because the VCF has duplicate IDs
VARID_FMT='@:#:$r:$a'

# --- 1. Full marker matrix for RKHS (no pruning) ----------------------------

plink2 --vcf "$VCF" \
  --export A \
  --out auspak_for_rkhs

# --- 2. LD pruning for BayesC ------------------------------------------------
# Window 50 SNPs, step 5, r2 threshold 0.5
# 1,267,773 / 1,824,377 variants removed -> ~557k markers retained

# Step 2a: identify markers to keep
plink2 --vcf "$VCF" \
  --set-all-var-ids "$VARID_FMT" \
  --indep-pairwise 50 5 0.5 \
  --out auspak_pruned_snps5

# Step 2b: export pruned marker matrix as mean-imputed .raw
plink2 --vcf "$VCF" \
  --set-all-var-ids "$VARID_FMT" \
  --extract auspak_pruned_snps5.prune.in \
  --export A \
  --out pruned05_AUSPAK_for_bayesC

# --- 3. Test subset (1000 random markers from the pruned set) ----------------

plink2 --vcf "$VCF" \
  --set-all-var-ids "$VARID_FMT" \
  --extract auspak_pruned_snps5.prune.in \
  --thin-count 1000 \
  --export A \
  --out AUSPAK_test_subset_1k