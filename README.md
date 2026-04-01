# Genomic Prediction in Quinoa: Comparing Statistical and Machine Learning Approaches

This repository contains code for multi-environment genomic prediction in quinoa (*Chenopodium quinoa*), comparing three statistical genetics methods (GBLUP, BayesC, RKHS) and a machine learning approach (LightGBM) across four cross-validation schemes.

## Overview

Genomic prediction was evaluated for 7 agronomic traits in 551 quinoa accessions grown across 6 environments (3 years each in Australia and Pakistan). Four models spanning different modelling paradigms were compared under four CV schemes that test different prediction scenarios relevant to breeding programs.

**Models:**

| Model | Package | Approach |
|-------|---------|----------|
| **GBLUP** | ASReml-R | Mixed model with genomic relationship matrix (VanRaden); location:year as random effect |
| **BayesC** | BGLR | Bayesian SNP regression with spike-and-slab prior; explicit marker effects |
| **RKHS** | BGLR | Bayesian kernel regression with Gaussian kernels at 3 bandwidths (kernel averaging) |
| **LightGBM** | LightGBM | Gradient boosting on principal components or kinship features; per-trait hyperparameter tuning |

**Cross-validation schemes:**

| Scheme | What it tests | Design |
|--------|---------------|--------|
| **CV1** | Predicting **new genotypes** in known environments | 5-fold GroupKFold on genotypes, 15 iterations |
| **CV2** | **Sparse testing** — random cells masked | 5-fold stratified by location-year, 15 iterations |
| **CV0** | **Leave-one-location-year-out** | Deterministic, entire environments held out |
| **CrossLoc** | **Cross-location transfer** | Train on one location, predict the other |

All models use the same fold assignments, seeds, preprocessing (z-score standardisation by location-year), and evaluation metrics for direct comparability.

## Repository Structure

```
.
├── README.md
├── preprocessing/
│   ├── SNPfiltering_for_AUSPAK_samples.sh   # VCF quality filtering (bcftools)
│   ├── prepare_marker_matrices.sh           # LD pruning + PLINK export
│   ├── Kinship_matrix.R                     # VanRaden kinship matrix (AGHmatrix)
│   ├── PCA.R                               # Principal components (SNPRelate)
│   └── BLUEs.ipynb                          # Phenotype processing (BLUEs)
│
├── data/
│   ├── AUSPAK_phenotypes_means_BLUEs.csv              # Phenotypes (BLUEs)
│   ├── AUSPAK_PCs_all.csv                             # 551 principal components
│   ├── kinship_matrix_VanRaden_auspak_maxmissing20.csv    # Kinship matrix (CSV)
│   ├── kinship_matrix_VanRaden_auspak_maxmissing20.RData  # Kinship matrix (RData)
│   ├── pruned05_AUSPAK_for_bayesC.raw                 # LD-pruned markers for BayesC (PLINK .raw)
│   ├── auspak_for_rkhs.raw                            # Full markers for RKHS (PLINK .raw)
│   └── AUSPAK_test_subset_1k.raw                      # 1000-marker test subset (551 genotypes)
│
├── models/
│   ├── GBLUP/                   # ASReml-R GBLUP pipeline
│   │   ├── G_matrix_GBLUP.R         # G matrix preparation (bend + invert)
│   │   ├── GBLUP_utils.R            # Shared utilities
│   │   └── GBLUP.R                  # All CV schemes in one script
│   │
│   ├── BayesC/                  # BGLR BayesC pipeline
│   │   ├── BayesC_utils.R           # Shared utilities
│   │   ├── BayesC_CV1_single_iter.R # CV1 (one iteration per SLURM job)
│   │   ├── BayesC_CV2_single_iter.R # CV2 (one iteration per SLURM job)
│   │   ├── BayesC_CV0.R             # CV0 (deterministic)
│   │   ├── BayesC_CrossLoc.R        # Cross-location prediction
│   │   ├── BayesC_full_model.R      # Full model diagnostics + marker effects
│   │   ├── launch_BayesC_jobs.sh    # SLURM job generation/submission
│   │   └── aggregate_BayesC_results.R
│   │
│   ├── RKHS/                    # BGLR RKHS kernel averaging pipeline
│   │   ├── RKHS_utils.R             # Shared utilities + kernel computation
│   │   ├── RKHS_CV1.R               # CV1 (all iterations per SLURM job)
│   │   ├── RKHS_CV2.R               # CV2 (all iterations per SLURM job)
│   │   ├── RKHS_CV0_CrossLoc.R      # CV0 + CrossLoc combined
│   │   ├── launch_RKHS_jobs.sh      # SLURM job generation/submission
│   │   └── aggregate_RKHS_results.R
│   │
│   ├── LightGBM/               # LightGBM gradient boosting pipeline
│   │   ├── Prepare_input_data.py    # Data preparation (3 pickle variants)
│   │   ├── tune_LightGBM.py         # Hyperparameter tuning (RandomizedSearchCV)
│   │   ├── LightGBM_utils.py        # Shared utilities
│   │   ├── run_LightGBM.py          # All CV schemes in one script
│   │   └── environment_ml.yml       # Conda environment specification
│   │
│   └── tests/                  # Unit + integration tests (R and Python)
│       ├── test_cv_fold_assignment.R     # Cross-model fold assignment tests
│       ├── test_*.R                      # Per-model R tests (BayesC, RKHS, GBLUP)
│       └── test_*.py                     # LightGBM Python tests
│
└── results/
    ├── combine_all_models.R         # Merge per-model results into unified CSVs
    ├── visualize_results.py         # Generate comparison figures (PDF)
    ├── cv_results_all_models.csv    # Combined accuracy metrics across all models
    └── cv_summary_all_models.csv    # Summary statistics across all models
```

## Data

### Phenotypic Traits

All traits are Best Linear Unbiased Estimates (BLUEs):

| Code | Trait | Direction |
|------|-------|-----------|
| DTF | Days to flowering | Lower is better |
| DTH | Days to harvest maturity | Lower is better |
| PtHt | Plant height (cm) | Lower is better |
| PcleLng | Panicle length (cm) | Higher is better |
| SdLen | Seed length (mm) | Higher is better |
| TGW | Thousand grain weight (g) | Higher is better |
| SdW_z | Seed yield (z-transformed) | Higher is better |

### Environments

| Location | Years | Location-years |
|----------|-------|----------------|
| AUS (Kununurra, Australia) | 2017, 2018, 2019 | AUS_2017, AUS_2018, AUS_2019 |
| PAK (Faisalabad, Pakistan) | 2019-20, 2020-21, 2021-22 | PAK_2019, PAK_2020, PAK_2021 |

551 quinoa accessions total. Not all genotypes appear in all location-years; not all traits are observed for every genotype x location-year combination.

### Genotypic Data

- 18 nuclear chromosomes (9 homeologous pairs: Cq1A-Cq9A, Cq1B-Cq9B)
- 1,824,377 biallelic SNPs after quality filtering (MAF >= 0.01, <20% missing, depth 5-30x)
- LD pruning (PLINK2, r² < 0.5, window 50, step 5) retains ~557k markers for BayesC; RKHS and GBLUP use the full marker set
- Kinship matrix: VanRaden method 1 (AGHmatrix)
- PCA: 551 principal components (SNPRelate)

## Models

### GBLUP

Mixed model using ASReml-R with a pre-computed genomic relationship matrix:

- **Fixed**: location
- **Random**: `vm(sample.id, Ginv_sparse)` (genomic BLUPs) + `location:year`
- G matrix is bent for positive definiteness and inverted to sparse triplet format (ASRgenomics)
- CrossLoc variant uses intercept-only fixed effects (`trait ~ 1 + G`)
- Runs locally (no SLURM needed)

### BayesC

Bayesian SNP regression via BGLR with a spike-and-slab prior:

- **Fixed**: location + year-within-location (as fixed effects, since BGLR lacks random effects)
- **Marker effects**: BayesC prior (some markers shrunk to zero)
- MCMC: 15,000 iterations, 5,000 burn-in, thinning every 5
- Full model mode provides variance components, Manhattan plots, and MCMC diagnostics
- Distributed via SLURM (one job per trait per iteration for CV1/CV2)

### RKHS

Reproducing Kernel Hilbert Space regression via BGLR with Gaussian kernel averaging:

- **Fixed**: location + year-within-location (same as BayesC)
- **Kernels**: 3 Gaussian kernels at different bandwidths (`h = 1/5, 1, 5` x `1/median(D)`)
- BGLR estimates variance components per kernel (Bayesian kernel averaging)
- Non-parametric: predicts through genomic similarity, not individual marker effects
- Distributed via SLURM (one job per trait)

### LightGBM

Gradient boosting with per-trait hyperparameter tuning:

- **Features**: genetic features (PCs or kinship columns) + one-hot encoded location and location-year
- Three input variants: all 551 PCs, 25 PCs, or kinship matrix columns
- Hyperparameters tuned via `RandomizedSearchCV` (50 combinations, 5-fold GroupKFold)
- CrossLoc uses genetic features only (no location encoding)
- Runs locally (no SLURM needed)

## Preprocessing

All preprocessing scripts are in `preprocessing/`. These must be run before the prediction models.

### 1. SNP Filtering

Quality filtering of the raw VCF (bcftools):

- Nuclear chromosomes only (excludes organellar genomes)
- Biallelic SNPs only
- MAF >= 0.01
- Missing data < 20% per variant
- Mean depth 5-30x (excludes low coverage and likely paralogs)

**Output:** `quinoa_551accessions_genomic_prediction.vcf` (1,824,377 SNPs, 551 accessions)

### 2. Marker Matrix Export and LD Pruning

`prepare_marker_matrices.sh` (PLINK2) exports marker data in three forms:

| Output file | Description | Used by |
|-------------|-------------|---------|
| `auspak_for_rkhs.raw` | Full marker set (1,824,377 SNPs) | RKHS |
| `pruned05_AUSPAK_for_bayesC.raw` | LD-pruned (r² < 0.5, ~557k markers) | BayesC |
| `AUSPAK_test_subset_1k.raw` | 1,000 random markers from pruned set | Tests |

### 3. Kinship Matrix

`Kinship_matrix.R` computes the VanRaden method 1 additive relationship matrix from the full VCF using AGHmatrix. Used by GBLUP (via G-inverse) and LightGBM (as kinship features).

### 4. Principal Component Analysis

`PCA.R` computes all 551 principal components from the full marker set using SNPRelate. Used by LightGBM as an alternative feature set.

### 5. Phenotype Processing (BLUEs)

`BLUEs.ipynb` estimates Best Linear Unbiased Estimates for each trait using ASReml-R, accounting for year and replicate effects. Outputs `AUSPAK_phenotypes_means_BLUEs.csv`.

### 6. Z-score Standardisation (within models)

All models standardise each trait within each location-year before fitting:

```
z_ijl = (y_ijl - mean_l) / sd_l
```

This removes environmental mean differences so that prediction accuracy reflects the ability to rank genotypes within environments, not to predict absolute trait values.

## Evaluation Metrics

All predictions are evaluated **per location-year** using:

- **Pearson correlation** — linear predictive accuracy
- **Spearman rank correlation** — rank-based predictive accuracy
- **NDCG@10** — normalised discounted cumulative gain at top 10, reflecting selection ranking quality at ~20% selection intensity

Trait directionality is accounted for (e.g., for DTF, lower values are ranked higher in NDCG).

## Usage

### Preprocessing (run once)

```bash
cd preprocessing
bash SNPfiltering_for_AUSPAK_samples.sh    # Quality filter VCF
bash prepare_marker_matrices.sh            # LD pruning + PLINK export
Rscript Kinship_matrix.R                   # Kinship matrix
Rscript PCA.R                              # Principal components
# Run BLUEs.ipynb in R/Jupyter             # Phenotype BLUEs
```

### R models (GBLUP, BayesC, RKHS)

```bash
# GBLUP — run locally
cd models/GBLUP
Rscript G_matrix_GBLUP.R          # Prepare G-inverse (once)
Rscript GBLUP.R                   # Run all CV schemes

# BayesC — SLURM cluster
cd models/BayesC
bash launch_BayesC_jobs.sh submit             # Submit all jobs
Rscript aggregate_BayesC_results.R all        # Aggregate after completion

# RKHS — SLURM cluster
cd models/RKHS
bash launch_RKHS_jobs.sh submit               # Submit all jobs
Rscript aggregate_RKHS_results.R all          # Aggregate after completion
```

### LightGBM (Python)

```bash
cd models/LightGBM
conda env create --file environment_ml.yml --prefix ./env-ML
conda activate ./env-ML

python Prepare_input_data.py       # Prepare pickles (once)
python tune_LightGBM.py            # Tune hyperparameters (once per input type)
python run_LightGBM.py             # Run all CV schemes
```

### Combine results and visualise

```bash
cd results
Rscript combine_all_models.R       # Merge all model results
python visualize_results.py        # Generate comparison figures (PDF)
```

## Tests

Unit and integration tests for all four model pipelines live in `models/tests/`:

```bash
cd models/tests

# R tests (BayesC, RKHS, GBLUP)
Rscript -e 'testthat::test_dir(".")'

# Python tests (LightGBM)
conda activate ../LightGBM/env-ML
pytest . -v
```

Tests cover fold assignment consistency, evaluation metrics, data loading, per-location-year evaluation, and end-to-end CV pipeline integration.

## Requirements

### External tools

- **bcftools** (v1.10+) — VCF filtering
- **PLINK2** — LD pruning and marker matrix export

### R packages

- `asreml` (requires license) + `ASRgenomics` — GBLUP and BLUEs
- `BGLR` — BayesC, RKHS
- `AGHmatrix` — kinship matrix
- `SNPRelate` — PCA
- `dplyr`, `data.table`, `Matrix` — data handling
- `testthat` — testing

### Python packages

See `models/LightGBM/environment_ml.yml`:
- `lightgbm`, `scikit-learn`, `pandas`, `numpy`, `scipy`
- `matplotlib` — visualisation
- `pytest` — testing

## Citation

If you use this code or data, please cite: (Manuscript currently under review)



## Contact

clara.stanschewski@kaust.edu.sa

