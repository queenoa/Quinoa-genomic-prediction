# BayesC Genomic Prediction Pipeline

## Overview

BayesC genomic prediction pipeline for the AUSPAK quinoa dataset, using the [BGLR](https://cran.r-project.org/package=BGLR) R package. Evaluates prediction accuracy for 7 traits across multiple environments (location-years) using four cross-validation schemes, plus a full model for diagnostics and marker effect estimation.

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
| **CV1** | `BayesC_CV1_single_iter.R` | Predicting **new genotypes** in known environments. 5-fold CV on genotypes, repeated across iterations (seed = 1000+iter). |
| **CV2** | `BayesC_CV2_single_iter.R` | **Sparse testing** — random observation-level cells masked. 5-fold CV stratified by location-year (seed = 2000+iter). |
| **CV0** | `BayesC_CV0.R` | **Leave-one-location-year-out** — entire environments held out. Deterministic, no random folds. Incremental CSV saving. |
| **CrossLoc** | `BayesC_CrossLoc.R` | **Cross-location transfer** — train on one location, predict all others using estimated marker effects (GEBVs). Incremental CSV saving. |

## Full model

`BayesC_full_model.R` trains on all data (no CV masking) for one trait. Outputs per trait in `full_model_<trait>/`:

| File | Contents |
|------|----------|
| `variance_components_<trait>.csv` | varE, varB, probIn, h2_marker, DIC, MCMC settings |
| `marker_effects_<trait>.csv` | Per-SNP posterior mean effect, sorted by |effect| |
| `fixed_effects_<trait>.csv` | Location + year-within-location estimates |
| `fitted_values_<trait>.csv` | Observed vs fitted for all genotype × location-year |
| `mcmc_chains_<trait>.rds` | Raw MCMC chain samples for convergence diagnostics |
| `diagnostics_<trait>.pdf` | Trace + density plots (varE, varB, mu, probIn) + Manhattan plot |

## Model structure

Defined in `BayesC_utils.R`. Each BGLR fit includes:

- **Fixed effects**: location + year-within-location (approximates an ASReml `fixed = trait ~ location`, `random = ~ location:year` design, but all fixed since BGLR lacks random effects)
- **BayesC marker effects**: SNP-level regression with a spike-and-slab prior (allows some markers to have zero effect)
- **Preprocessing**: z-score standardisation of the trait within each location-year; mean-imputation of missing marker data

MCMC settings: 15,000 iterations, 5,000 burn-in, thinning every 5. Overridable via environment variables `BAYESC_NITER`, `BAYESC_BURNIN`, `BAYESC_THIN`.

## Evaluation metrics

Each CV script saves a companion predictions CSV (`predictions_<scheme>_<trait>_...BayesC.csv`) with one row per test accession × location-year, recording `sample.id`, `location_year`, `location`, `observed`, and `predicted` values. These allow post-hoc recomputation of any accuracy metric, examination of prediction bias, or identification of outlier genotypes.

Predictions are evaluated **per location-year** using:

- **Pearson correlation** — linear predictive accuracy
- **Spearman rank correlation** — rank-based predictive accuracy
- **NDCG@10** — normalized discounted cumulative gain at top 10, for selection ranking quality

## Files

| File | Purpose |
|------|---------|
| `BayesC_utils.R` | Shared utilities: data loading, z-score scaling, marker matrix construction, BGLR model fitting wrapper, evaluation metrics |
| `BayesC_CV1_single_iter.R` | Runs one iteration of CV1 (new genotypes). Called per trait per iteration. |
| `BayesC_CV2_single_iter.R` | Runs one iteration of CV2 (sparse testing). Called per trait per iteration. |
| `BayesC_CV0.R` | Runs CV0 for one trait (deterministic, incremental CSV saving). |
| `BayesC_CrossLoc.R` | Runs CrossLoc for one trait (deterministic, incremental CSV saving). |
| `BayesC_CV0_CrossLoc.R` | Legacy combined CV0 + CrossLoc script (kept for reference). |
| `BayesC_full_model.R` | Trains on all data — variance components, marker effects, MCMC diagnostics. |
| `launch_BayesC_jobs.sh` | Generates and optionally submits SLURM jobs. Supports iteration ranges and selective mode submission. |
| `aggregate_BayesC_results.R` | Collects per-iteration CSV outputs (accuracy + predictions) into combined results, summary statistics, and a combined predictions file after all jobs complete. |

## Usage

```bash
# Generate SLURM scripts only
bash launch_BayesC_jobs.sh

# Generate and submit all traits (all CV schemes)
bash launch_BayesC_jobs.sh submit

# Submit specific iteration range (CV1 + CV2 only)
bash launch_BayesC_jobs.sh submit 1 2
bash launch_BayesC_jobs.sh submit DTF_blue 3 5

# Submit specific CV schemes
bash launch_BayesC_jobs.sh submit cv0
bash launch_BayesC_jobs.sh submit DTF_blue crossloc

# Submit full model (diagnostics + marker effects)
bash launch_BayesC_jobs.sh submit fullmodel
bash launch_BayesC_jobs.sh submit DTF_blue fullmodel

# After all CV jobs finish, aggregate results
Rscript aggregate_BayesC_results.R all
```

## Input data

### Environments

- **Locations (2):** `AUS` (Kununurra, Australia), `PAK` (Faisalabad, Pakistan)
- **Location-years (6):** `AUS_2017`, `AUS_2018`, `AUS_2019`, `PAK_2019`, `PAK_2020`, `PAK_2021`
- Not all genotypes appear in all location-years; not all traits are observed for every genotype × location-year combination (NAs present)

### Raw data files

- **Phenotype file**: `AUSPAK_phenotypes_means_BLUEs.csv` — must contain columns `sample.id`, `location`, `year`, `location_year`, and the trait column
- **Marker file**: `pruned05_AUSPAK_for_bayesC.raw` — PLINK `.raw` format with IID column and SNP dosages (first 6 columns are PLINK metadata)
- **Test data**: `data/AUSPAK_test_subset_1k.raw` — 1000-marker subset (551 genotypes) for running tests

## Tests

Unit and integration tests in `tests/` using `testthat`. Run from the `tests/` directory:

```bash
cd tests
Rscript -e 'testthat::test_dir(".")'
```

Test files cover: `calculate_ndcg`, `evaluate_predictions`, `apply_location_year_scaling`, `build_fixed_design`, `build_obs_marker_matrix_and_groups`, `load_and_prepare_data`, `evaluate_per_location_year`, `aggregate_helpers`, and end-to-end integration of CV scripts.
