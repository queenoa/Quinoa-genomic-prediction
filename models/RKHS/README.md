# RKHS Multi-Kernel Genomic Prediction Pipeline

## Overview

Reproducing Kernel Hilbert Space (RKHS) genomic prediction pipeline for the AUSPAK quinoa dataset, using the [BGLR](https://cran.r-project.org/package=BGLR) R package. Evaluates prediction accuracy for 7 traits across multiple environments (location-years) using four cross-validation schemes.

Where BayesC models explicit SNP marker effects, RKHS captures genomic relationships through Gaussian kernels computed from the full marker matrix. The kernel approach is non-parametric: it does not estimate individual marker effects but instead predicts through similarity in the marker space. Both pipelines share the same fixed-effect structure, preprocessing, evaluation metrics, and CV designs.

## RKHS method

### Gaussian kernel construction

1. Load PLINK `.raw` marker matrix, mean-impute missing genotypes
2. Center and scale each SNP (zero mean, unit variance) so all markers contribute equally regardless of MAF
3. Compute squared Euclidean distance matrix: `D[i,j] = ||x_i - x_j||^2`
4. Normalise: `D <- D / mean(D)`
5. Compute three Gaussian kernels at different bandwidths: `K(x_i, x_j) = exp(-h * D[i,j])`
   - `h1 = 1/(5 * median(D))` — wide (smooth, captures broad similarity)
   - `h2 = 1/median(D)` — medium
   - `h3 = 5/median(D)` — narrow (captures fine-grained similarity)

Following Perez & de los Campos (2014, Box 11) and Cuevas et al. (2016). BGLR estimates the variance component for each kernel separately, effectively performing Bayesian model averaging over bandwidths (kernel averaging, KA).

### Kernel expansion to observation level

For multi-environment data in long format, observation `i` corresponds to genotype `g_i`. The observation-level kernel is simply: `K_obs[i,j] = K_geno[g_i, g_j]`. This is the kernel analogue of `build_obs_marker_matrix()` in BayesC, but much cheaper (indexing vs matrix multiplication).

Kernels are computed once from the full marker matrix and checkpointed to `RKHS_kernels.RData`.

### Key difference from BayesC: empty location-years

Unlike BayesC, empty location-years (100% missing for the target trait) are **not** dropped. The kernel is precomputed at the genotype level and expanded to observation level by simple indexing, so the computational cost of extra NA rows is negligible.

## Traits

All input traits are BLUEs (Best Linear Unbiased Estimates):

| Code | Trait |
|------|-------|
| `DTF_blue` | Days to flowering |
| `DTH_blue` | Days to maturity |
| `PtHt_blue` | Plant height |
| `PcleLng_blue` | Panicle length |
| `SdLen_blue` | Seed length |
| `TGW_blue` | Thousand grain weight |
| `SdW_z_blue` | Seed yield |

Traits where lower is better (DTF, DTM, plant height) have NDCG@10 sign-flipped accordingly.

## Cross-validation schemes

| Scheme | Script | What it tests |
|--------|--------|---------------|
| **CV1** | `RKHS_CV1.R` | Predicting **new genotypes** in known environments. 5-fold CV on genotypes, repeated across 15 iterations (seed = 1000+iter). |
| **CV2** | `RKHS_CV2.R` | **Sparse testing** — random observation-level cells masked. 5-fold CV stratified by location-year, 15 iterations (seed = 2000+iter). |
| **CV0** | `RKHS_CV0_CrossLoc.R` | **Leave-one-location-year-out** — entire environments held out. Deterministic, no random folds. |
| **CrossLoc** | `RKHS_CV0_CrossLoc.R` | **Cross-location transfer** — train on one location, predict all others via kernel similarity. |

All four schemes match the BayesC pipeline exactly in design (fold assignment, seed scheme, stratification).

## Model structure

Defined in `RKHS_utils.R`. Each BGLR fit includes:

- **Fixed effects**: location + year-within-location, via `build_fixed_design()` (identical to BayesC; approximates ASReml `fixed = trait ~ location`, `random = ~ location:year`, but all fixed)
- **3 RKHS kernels**: Gaussian kernels at narrow, medium, and wide bandwidths — BGLR estimates variance components for each, performing Bayesian kernel averaging
- **Single residual variance**: BGLR does not support heterogeneous residual variances (`groups` argument) with RKHS model type — unlike BayesC, RKHS uses one residual variance across all environments
- **Preprocessing**: z-score standardisation of the trait within each location-year; mean-imputation of missing marker data

MCMC settings: 15,000 iterations, 5,000 burn-in, thinning every 5. Overridable via environment variables `RKHS_NITER`, `RKHS_BURNIN`, `RKHS_THIN`.

### CrossLoc model variant

For cross-location prediction, the model is simpler:
- **No location fixed effect** (target location cannot be estimated)
- **Year fixed effect within training location** (if >1 year), with target rows set to the reference level
- **3 RKHS kernels** built from the combined (train + target) genotype set
- Target genotypes have `y = NA`, so BGLR predicts them from kernel similarity to training genotypes

This is the kernel analogue of the BayesC CrossLoc approach (which extracts marker effects and computes GEBVs). Here, prediction flows through the kernel rather than explicit marker effects.

## Evaluation metrics

Each CV script saves a companion predictions CSV (`predictions_<scheme>_<trait>_RKHS.csv`) with one row per test observation, recording `sample.id`, `location_year`, `location`, `observed`, and `predicted` values. These allow post-hoc recomputation of any accuracy metric.

Predictions are evaluated **per location-year** using:

- **Pearson correlation** — linear predictive accuracy
- **NDCG@10** — normalized discounted cumulative gain at top 10, for selection ranking quality

## Files

| File | Purpose |
|------|---------|
| `RKHS_utils.R` | Shared utilities: data loading, kernel computation/checkpointing, z-score scaling, observation-level kernel expansion, fixed-effect design matrix, BGLR RKHS fitting wrapper, evaluation metrics |
| `RKHS_CV1.R` | Runs all 15 iterations of CV1 (new genotypes). One SLURM job per trait. |
| `RKHS_CV2.R` | Runs all 15 iterations of CV2 (sparse testing, 5-fold stratified). One SLURM job per trait. |
| `RKHS_CV0_CrossLoc.R` | Runs both CV0 and CrossLoc for one trait (deterministic, incremental CSV saving). |
| `launch_RKHS_jobs.sh` | Generates and optionally submits SLURM jobs (3 per trait). |
| `aggregate_RKHS_results.R` | Collects per-scheme CSV outputs (accuracy + predictions) into combined results and summary statistics. |

## Usage

```bash
# Generate SLURM scripts only
bash launch_RKHS_jobs.sh

# Generate and submit all traits
bash launch_RKHS_jobs.sh submit

# Submit one trait only
bash launch_RKHS_jobs.sh submit DTF_blue

# After all jobs finish, aggregate results
Rscript aggregate_RKHS_results.R all
Rscript aggregate_RKHS_results.R DTF_blue
```

## Input data

### Environments

- **Locations (2):** `AUS` (Kununurra, Australia), `PAK` (Faisalabad, Pakistan)
- **Location-years (6):** `AUS_2017`, `AUS_2018`, `AUS_2019`, `PAK_2019`, `PAK_2020`, `PAK_2021`
- Not all genotypes appear in all location-years; not all traits are observed for every genotype × location-year combination (NAs present)

### Raw data files

- **Phenotype file**: `AUSPAK_phenotypes_means_BLUEs.csv` — must contain columns `sample.id`, `location`, `year`, `location_year`, and the trait column
- **Marker file**: `auspak_for_rkhs.raw` — PLINK `.raw` format with IID column and SNP dosages (first 6 columns are PLINK metadata)
- **Kernel checkpoint**: `RKHS_kernels.RData` — computed once from the marker file and reused across all traits and CV schemes

## Output files

Per-trait, per-scheme:
- `cv_results_<scheme>_<trait>_RKHS.csv` — accuracy metrics per location-year
- `predictions_<scheme>_<trait>_RKHS.csv` — individual-level predictions

After aggregation:
- `cv_results_<trait>_RKHS_all_schemes.csv` — combined row-level results
- `cv_summary_<trait>_RKHS_all_schemes.csv` — summary statistics (mean, SD, min, max per scheme per location)
- `predictions_<trait>_RKHS_all_schemes.csv` — combined individual predictions

## Functions shared with BayesC

The following functions in `RKHS_utils.R` are identical to their counterparts in `BayesC_utils.R`:

- `VALID_TRAITS`, `LOWER_IS_BETTER_TRAITS` — trait constants
- `apply_location_year_scaling()` — z-score standardisation by location-year
- `build_fixed_design()` — fixed-effect design matrix (location + year-within-location)
- `build_groups()` — environment group vector (defined but not used in RKHS fitting; BGLR does not support `groups` with RKHS model type)
- `calculate_ndcg()` — NDCG@k metric
- `evaluate_predictions()` — Pearson, Spearman, NDCG@10 evaluation
- `evaluate_per_location_year()` — per-location-year evaluation loop


