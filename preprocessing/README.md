# Preprocessing

Builds the marker, kinship, and phenotype inputs consumed by the genomic
prediction models in `../models/`. The top-level `../README.md` gives the
project-wide context; this README documents the scripts in this directory.

## Pipeline

Run in this order — each step depends on outputs from the previous one(s):

1. `SNPfiltering_for_AUSPAK_samples.sh` — quality-filter the raw VCF
2. `prepare_marker_matrices.sh` — export PLINK `.raw` marker matrices
3. `Kinship_matrix.R` — compute the VanRaden kinship matrix
4. `build_GP_phenotype_input.R` — build per-location-year phenotype means

The `BLUEs/` subfolder is a separate, optional pipeline used only for
phenotype visualisation (see below).

## Scripts

### `SNPfiltering_for_AUSPAK_samples.sh`

Quality filtering of the raw VCF using `bcftools` (SLURM batch script).
Three steps in one job:

1. Subset the master VCF to the 551 AUS/PAK accessions
   (`auspak_sample_list.txt`).
2. Apply variant QC filters:
   - Nuclear chromosomes only (`Cq1A–Cq9B`; excludes organellar genomes)
   - Biallelic SNPs only (`-m2 -M2 -v snps`)
   - `F_MISSING < 0.2` (≤ 20% missing per variant)
   - `MAF ≥ 0.01`
   - Mean depth 5–30× (excludes low coverage and likely paralogs)
3. Strip everything except the `GT` field and drop GATK command lines from
   the header.

**Output:** `quinoa_551accessions_genomic_prediction.vcf` (1,824,377 SNPs,
551 accessions).

### `prepare_marker_matrices.sh`

PLINK2 export of the filtered VCF into three `.raw` (mean-imputed)
matrices. Variant IDs are reassigned to `chr:pos:ref:alt` because the
source VCF has duplicate IDs.

| Output | Description | Used by |
|--------|-------------|---------|
| `auspak_for_rkhs.raw` | Full marker set (1,824,377 SNPs) | RKHS |
| `pruned05_AUSPAK_for_bayesC.raw` | LD-pruned (`--indep-pairwise 50 5 0.5`, ~557k SNPs) | BayesC |
| `AUSPAK_test_subset_1k.raw` | 1,000 random markers from the pruned set | `models/tests/` |

### `Kinship_matrix.R`

Computes the VanRaden method 1 additive relationship matrix from the
filtered VCF using `AGHmatrix::Gmatrix`. Genotype calls are recoded
`0/0 → 0`, `0/1 → 1`, `1/1 → 2` (both phased and unphased separators
handled). The missing-data threshold is set to 20% to match the variant
filter applied upstream.

**Output:** `kinship_matrix_VanRaden_auspak_maxmissing20.RData`. Consumed
by GBLUP (as the G matrix) and LightGBM (as kinship features).

### `build_GP_phenotype_input.R`

Builds the phenotype input file for every genomic prediction model:

- **Output:** `../data/AUSPAK_phenotypes_GP_input.csv`
- **Format:** wide — one row per (location, year, accession), one column
  per trait (`DTF`, `DTH`, `PtHt`, `PcleLng`, `SdLen`, `TGW`, `SdW_z`).
- **Source:** raw per-accession means from `../data/AUS_phenotypes_raw.csv`
  and `../data/PAK_phenotypes_raw.csv`.
- **Filter:** only sequenced accessions (non-NA `SampleName`) are kept.

Raw means are used uniformly across all six location-years. Five of the
six trials are unreplicated (mean = single observation). PAK 2021-22 is
the only replicated trial (2 reps), but without spatial correction or
additional covariates the BLUEs are equal to raw means for the fully
replicated accessions and the adjustment for the few partially replicated
ones is negligible — so means are used everywhere to keep the GP input
uniform.

## `BLUEs/` (visualisation only)

`BLUEs/` fits per-location BLUEs (one value per accession per location per
trait) for use in phenotype-correlation plots and descriptive summaries.
These BLUEs are **not** consumed by the genomic prediction models. See
`BLUEs/README.md` for details.

## Usage

```bash
# From this directory:
bash SNPfiltering_for_AUSPAK_samples.sh   # SLURM script — submit with sbatch
bash prepare_marker_matrices.sh           # local
Rscript Kinship_matrix.R                  # local
Rscript build_GP_phenotype_input.R        # local

# Optional — for phenotype-correlation plots only:
cd BLUEs
Rscript BLUEs_AUS.R
Rscript BLUEs_PAK.R
Rscript merge_BLUEs_for_phenotype_plots.R
```

## Requirements

- `bcftools` (v1.10+)
- `plink2`
- R packages: `data.table`, `AGHmatrix`, `dplyr` (and `asreml` +
  `ASRgenomics` for the optional `BLUEs/` pipeline)
