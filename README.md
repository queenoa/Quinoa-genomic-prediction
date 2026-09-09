# Genomic Prediction in Quinoa: Comparing Statistical and Machine Learning Approaches

This repository contains code for multi-environment genomic prediction in quinoa (*Chenopodium quinoa*), comparing three statistical genetics methods (GBLUP, BayesC, RKHS) and a machine-learning approach (LightGBM) across four cross-validation schemes.

## Overview

Genomic prediction was evaluated for 7 agronomic traits in 551 quinoa accessions grown across 6 environments (3 years each in Australia and Pakistan). Four models spanning different modelling paradigms were compared under four CV schemes that test different prediction scenarios relevant to breeding programs.

**Models:**

| Model | Package | Approach |
|-------|---------|----------|
| **GBLUP** | ASReml-R | Mixed model with VanRaden genomic relationship matrix; heterogeneous residual variance per location; location-specific year effects. Also fit per-location (within each location independently) to benchmark and estimate location-specific variance components. |
| **BayesC** | BGLR | Bayesian SNP regression with spike-and-slab prior; explicit marker effects |
| **RKHS** | BGLR | Bayesian kernel regression with Gaussian kernels at 3 bandwidths (kernel averaging) |
| **LightGBM** | LightGBM | Gradient boosting on kinship features + one-hot-encoded environments; per-trait hyperparameter tuning |

**Cross-validation schemes:**

| Scheme | What it tests | Design |
|--------|---------------|--------|
| **CV1** | Predicting **new genotypes** in known environments | 5-fold on genotypes, 15 iterations |
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
│   ├── build_GP_phenotype_input.R           # Builds AUSPAK_phenotypes_GP_input.csv (per-location means)
│   └── BLUEs/                               # Per-location BLUEs for visualisation only (asreml-R)
│
├── data/
│   ├── AUSPAK_phenotypes_GP_input.csv             # Phenotypes for genomic prediction (per-location means)
│   ├── AUSPAK_BLUEs_for_phenotype_plots.csv       # Per-location BLUEs for phenotype-correlation plots
│   ├── kinship_matrix_VanRaden_auspak_maxmissing20.csv    # Kinship matrix (CSV)
│   ├── kinship_matrix_VanRaden_auspak_maxmissing20.RData  # Kinship matrix (RData)
│   ├── pruned05_AUSPAK_for_bayesC.raw       # LD-pruned markers for BayesC (PLINK .raw)
│   ├── auspak_for_rkhs.raw                  # Full markers for RKHS (PLINK .raw)
│   └── AUSPAK_test_subset_1k.raw            # 1000-marker test subset (551 genotypes)
│
├── models/
│   ├── GBLUP/        # ASReml-R GBLUP (global + per-location)
│   ├── BayesC/       # BGLR BayesC
│   ├── RKHS/         # BGLR RKHS kernel averaging
│   ├── LightGBM/     # LightGBM gradient boosting (kinship features)
│   └── tests/        # Unit + integration tests (R + Python)
│
└── visualizations_results/
    ├── combine_all_models.R                       # Merge per-model results into unified CSVs
    ├── Manuscript_Figures.ipynb                   # Notebook producing every manuscript figure
    ├── cv_results_all_models.csv                  # Combined accuracy metrics across all models (global CV1/CV2/CV0/CrossLoc)
    ├── cv_results_cv1_per-location-model.csv     # Per-location GBLUP CV1 accuracies
    └── Figure1.png … Figure7.png, FigureS1.png    # Rendered manuscript figures
```

Each model directory has a `README.md` with the implementation details, file inventory, and usage commands specific to that pipeline.

## Data

### Phenotypic Traits

Input data for the genomic prediction models consist of location-level phenotypic means (`AUSPAK_phenotypes_GP_input.csv`), not BLUEs. Trial-level means were used instead of BLUEs to retain environment-specific information for cross-validation. Owing to the largely unreplicated trial structure — with formal replication present only within one of the six location-year environments — insufficient information was available to estimate reliable within-environment BLUEs. Even in the trial containing two replicates, BLUEs provided a limited advantage over unadjusted means in the absence of additional covariates.

Per-location BLUEs (`AUSPAK_BLUEs_for_phenotype_plots.csv`), produced by the scripts in `preprocessing/BLUEs/`, are kept for visualisation only (e.g. trait correlation plots).

| Code | Trait | Direction |
|------|-------|-----------|
| `DTF` | Days to flowering | Lower is better |
| `DTH` | Days to harvest maturity | Lower is better |
| `PtHt` | Plant height (cm) | Lower is better |
| `PcleLng` | Panicle length (cm) | Higher is better |
| `SdLen` | Seed length (mm) | Higher is better |
| `TGW` | Thousand grain weight (g) | Higher is better |
| `SdW_z` | Seed yield (z-transformed) | Higher is better |


### Environments

| Location | Years | Location-years |
|----------|-------|----------------|
| AUS (Kununurra, Australia) | 2017, 2018, 2019 | AUS_2017, AUS_2018, AUS_2019 |
| PAK (Faisalabad, Pakistan) | 2019-20, 2020-21, 2021-22 | PAK_2019, PAK_2020, PAK_2021 |

551 quinoa accessions total. Not all genotypes appear in all location-years; not all traits are observed for every genotype × location-year combination.

### Genotypic Data

- 18 nuclear chromosomes (9 homeologous pairs: Cq1A–Cq9A, Cq1B–Cq9B)
- 1,824,377 biallelic SNPs after quality filtering (MAF ≥ 0.01, < 20% missing, depth 5–30×)
- LD pruning (PLINK2, r² < 0.5, window 50, step 5) retains ~557k markers — used **only for BayesC**; RKHS uses the full 1,824,377-SNP quality-filtered set
- Kinship matrix: VanRaden method 1 (AGHmatrix), built from the **full quality-filtered SNP set** (not the LD-pruned subset) — used by GBLUP (as G-inverse) and LightGBM (as kinship features)

## Models

### GBLUP

Mixed model using ASReml-R with a pre-computed genomic relationship matrix. Two variants:

- **Global model** (pooled across both locations) — `fixed = ~ location`, `random = vm(sample.id, G) + at(location):year` (location-specific year variance components), `residual = dsum(~ units | location)` (heterogeneous residual variance per location). Used for CV1, CV2, CV0. The CrossLoc variant uses an intercept-only fixed effect with `~ units` residual since `dsum`/`at(location)` collapse to homogeneous within a single location.
- **Per-location model** — fit independently within each location: `fixed = ~ 1`, `random = vm(sample.id, G) + year`, `residual = ~ units`. Used to benchmark within-location vs. global prediction accuracy and to estimate location-specific genetic variance components. Heterogeneous residual variances across years are not modelled at the per-location scale due to convergence failures (2–3 years per location, no within-year replication, year already random).

Runs locally (no SLURM). See [models/GBLUP/README.md](models/GBLUP/README.md).

### BayesC

Bayesian SNP regression via BGLR with a spike-and-slab prior:

- **Fixed**: location + year-within-location (BGLR lacks random effects, so these are fitted as fixed)
- **Marker effects**: BayesC prior (some markers shrunk to zero)
- MCMC: 15,000 iterations, 5,000 burn-in, thinning every 5
- Distributed via SLURM (one job per trait per iteration for CV1/CV2)

See [models/BayesC/README.md](models/BayesC/README.md).

### RKHS

Reproducing Kernel Hilbert Space regression via BGLR with Gaussian kernel averaging:

- **Fixed**: location + year-within-location (same as BayesC)
- **Kernels**: 3 Gaussian kernels at different bandwidths (`h = 1/5, 1, 5` × `1/median(D)`)
- BGLR estimates variance components per kernel (Bayesian kernel averaging)
- Non-parametric: predicts through genomic similarity, not individual marker effects
- Distributed via SLURM (one job per trait)

See [models/RKHS/README.md](models/RKHS/README.md).

### LightGBM

Gradient boosting on the kinship matrix:

- **Features**: kinship matrix columns (`K_*`) + one-hot-encoded `location` and `location_year` (CrossLoc uses kinship only)
- **Hyperparameters**: tuned per trait via `RandomizedSearchCV` (50 combinations, 5-fold GroupKFold)
- Runs locally (no SLURM)

See [models/LightGBM/README.md](models/LightGBM/README.md).

## Preprocessing

All preprocessing scripts are in `preprocessing/`. Run before the prediction models. See `preprocessing/README.md` for full details.

### 1. SNP filtering

`SNPfiltering_for_AUSPAK_samples.sh` (bcftools, SLURM script) subsets the master VCF to the 551 AUS/PAK accessions and applies variant QC:

- Nuclear chromosomes only (`Cq1A–Cq9B`; excludes organellar genomes)
- Biallelic SNPs only
- MAF ≥ 0.01
- Missing data < 20% per variant
- Mean depth 5–30× (excludes low coverage and likely paralogs)

**Output:** `quinoa_551accessions_genomic_prediction.vcf` (1,824,377 SNPs, 551 accessions)

### 2. Marker matrix export and LD pruning

`prepare_marker_matrices.sh` (PLINK2) exports the filtered VCF as three mean-imputed `.raw` matrices. Variant IDs are reassigned to `chr:pos:ref:alt` to ensure no duplicate IDs. **LD pruning is applied only to the BayesC input**; RKHS and the kinship matrix use the full quality-filtered SNP set.

| Output file | Description | Used by |
|-------------|-------------|---------|
| `auspak_for_rkhs.raw` | Full quality-filtered marker set (1,824,377 SNPs) | RKHS |
| `pruned05_AUSPAK_for_bayesC.raw` | LD-pruned subset of the full set (`--indep-pairwise 50 5 0.5`, ~557k markers) | BayesC |
| `AUSPAK_test_subset_1k.raw` | 1,000 random markers from the pruned set | Tests |

### 3. Kinship matrix

`Kinship_matrix.R` computes the VanRaden method 1 additive relationship matrix from the **full quality-filtered VCF** (not the LD-pruned subset) using `AGHmatrix::Gmatrix`. Used by GBLUP (as the G matrix) and LightGBM (as kinship features).

**Output:** `kinship_matrix_VanRaden_auspak_maxmissing20.RData`

### 4. Phenotype input

`build_GP_phenotype_input.R` produces `AUSPAK_phenotypes_GP_input.csv` — per-location per-accession means across all trials, in wide format (one row per location-year-accession, one column per trait). This is the input file consumed by every genomic prediction model. Raw means are used uniformly across all six location-years: five trials are unreplicated and the only replicated trial (PAK 2021-22) yields BLUEs essentially equal to means in the absence of spatial covariates.

### 5. Per-location BLUEs (visualisation only)

`preprocessing/BLUEs/` fits per-location BLUEs (asreml-R) for phenotype-correlation plots and descriptive summaries — **not** consumed by any GP model. `BLUEs_AUS.R` and `BLUEs_PAK.R` fit per-location-year mixed models; `merge_BLUEs_for_phenotype_plots.R` merges them into `AUSPAK_BLUEs_for_phenotype_plots.csv`.

### 6. Z-score standardisation (within models)

All models standardise each trait within each location-year before fitting:

```
z_ijl = (y_ijl - mean_l) / sd_l
```

This removes environmental mean differences so prediction accuracy reflects the ability to rank genotypes within environments, not to predict absolute trait values.

## Evaluation Metrics

All predictions are evaluated **per location-year** using:

- **Pearson correlation** — linear predictive accuracy
- **Spearman rank correlation** — rank-based predictive accuracy
- **NDCG@10** — normalised discounted cumulative gain at top 10, reflecting selection ranking quality at ~20% selection intensity

Trait directionality is accounted for (e.g., for DTF lower values are ranked higher in NDCG).

## Usage

### Preprocessing (run once)

```bash
cd preprocessing
bash SNPfiltering_for_AUSPAK_samples.sh    # Quality filter VCF (SLURM)
bash prepare_marker_matrices.sh            # LD pruning + PLINK export
Rscript Kinship_matrix.R                   # Kinship matrix
Rscript build_GP_phenotype_input.R         # Per-location means → AUSPAK_phenotypes_GP_input.csv

# (optional) Per-location BLUEs for visualisation only
cd BLUEs
Rscript BLUEs_AUS.R
Rscript BLUEs_PAK.R
Rscript merge_BLUEs_for_phenotype_plots.R
```

### R models (GBLUP, BayesC, RKHS)

```bash
# GBLUP — run locally
cd models/GBLUP
Rscript G_matrix_GBLUP.R              # Prepare G-inverse (once)
Rscript GBLUP.R                       # Global CV pipeline (all four schemes)
Rscript GBLUP_per_location_CV1.R      # Per-location CV1 + variance components

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

python Prepare_input_data.py       # Build model_inputs/model_input.pkl (once)
python tune_LightGBM.py            # Tune hyperparameters (once)
python run_LightGBM.py             # Run all CV schemes
```

### Combine results and visualise

```bash
cd visualizations_results
Rscript combine_all_models.R       # Merge all model results → cv_results_all_models.csv
jupyter notebook Manuscript_Figures.ipynb   # Reproduce every manuscript figure (Figure1–7, FigureS1)
```

All code used to generate the manuscript figures lives in `Manuscript_Figures.ipynb`; the rendered figures (`Figure1.png` … `FigureS1.png`) are committed alongside it.

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

Tests cover fold assignment consistency, evaluation metrics, data loading, per-location-year evaluation, and end-to-end CV pipeline integration. See [models/tests/README.md](models/tests/README.md).

## Requirements

### External tools

- **bcftools** (v1.10+) — VCF filtering
- **PLINK2** — LD pruning and marker matrix export

### R packages

- `asreml` (requires license) + `ASRgenomics` — GBLUP (and the optional `preprocessing/BLUEs/` pipeline)
- `BGLR` — BayesC, RKHS
- `AGHmatrix` — kinship matrix
- `dplyr`, `data.table`, `Matrix` — data handling
- `testthat` — testing

### Python packages

See `models/LightGBM/environment_ml.yml`:
- `lightgbm`, `scikit-learn`, `pandas`, `numpy`, `scipy`
- `matplotlib` — visualisation
- `pytest` — testing

## Citation

If you use this code or data, please cite: 
[Stanschewski, C. S., Warmington, M., Afzal, I., Rey, E., Fiene, G., Craine, E., ... & Poland, J. (2026). Genomic prediction in quinoa across contrasting environments using statistical and machine learning models. The Plant Genome, 19(3), e70277.](https://doi.org/10.1002/tpg2.70277)

## Contact

clara.stanschewski@kaust.edu.sa
