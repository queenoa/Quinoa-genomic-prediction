# Test Suite for AUSPAK Genomic Selection Models

## Overview

Shared test directory covering all four model pipelines: **BayesC**, **RKHS**, **GBLUP** (R, testthat), and **LightGBM** (Python, pytest). Tests range from lightweight unit tests of individual utility functions to end-to-end integration tests that fit models on a 1000-marker test subset.

LightGBM tests live separately under `models/LightGBM/tests/` because they depend on the Python environment and import directly from the LightGBM scripts.

## Running tests

### R tests (BayesC, RKHS, GBLUP)

```bash
cd tests

# Run all R tests
Rscript -e 'testthat::test_dir(".")'

# Run a single test file
Rscript -e 'testthat::test_file("test_cv_fold_assignment.R")'

# Run only tests matching a pattern
Rscript -e 'testthat::test_dir(".", filter = "gblup")'
```

**Dependencies:** `testthat`, `dplyr`, `data.table`, `BGLR` (BayesC/RKHS integration tests), `asreml` + `ASRgenomics` (GBLUP integration tests).

**Test data:** `data/AUSPAK_test_subset_1k.raw` (1000-marker subset, 551 genotypes) and `data/AUSPAK_phenotypes_means_BLUEs.csv`.

### Python tests (LightGBM)

```bash
cd models/LightGBM
conda activate ./env-ML
pytest tests/ -v
```

**Dependencies:** `pytest`, `scikit-learn`, `lightgbm`, `pandas`, `numpy`.

## Test speed tiers

| Tier | Files | Time | What runs |
|------|-------|------|-----------|
| **Instant** (<1s) | `test_cv_fold_assignment.R`, `test_calculate_ndcg.R`, `test_evaluate_predictions.R`, `test_apply_location_year_scaling.R`, `test_build_fixed_design.R`, `test_build_obs_marker_matrix_and_groups.R` | <1s each | Pure R logic, no model fitting |
| **Fast** (1-10s) | `test_evaluate_per_location_year.R`, `test_load_and_prepare_data.R`, `test_aggregate_helpers.R`, `test_gblup_shared_functions.R`, `test_gblup_genotype_alignment.R`, `test_gblup_evaluate_per_location_year.R`, `test_rkhs_shared_functions.R`, `test_rkhs_build_obs_kernel.R`, `test_rkhs_aggregate_helpers.R`, `test_rkhs_kernel_computation.R` | 1-10s each | Data loading, matrix ops, no model fitting |
| **Slow** (1-10 min) | `test_integration_cv_scripts.R`, `test_rkhs_integration_cv_scripts.R`, `test_gblup_integration_cv_functions.R`, `test_rkhs_fit_and_predict.R`, `test_gblup_fit_and_predict.R`, `test_rkhs_load_and_prepare_data.R` | 1-10 min each | BGLR/ASReml model fitting with minimal MCMC |

To run only the fast tests (no model fitting):

```bash
Rscript -e 'testthat::test_dir(".", filter = "^(?!.*integration|.*fit_and_predict|.*rkhs_load)")'
```

## Test files by model

### Cross-model (shared logic)

| File | Tests | What it covers |
|------|-------|----------------|
| `test_cv_fold_assignment.R` | 11 | CV1 and CV2 fold creation logic shared by all R models: shuffling across iterations, reproducibility, balance, stratification, unobserved row handling |
| `test_calculate_ndcg.R` | 10 | NDCG@k: perfect/worst rankings, lower-is-better, ties, negatives, edge cases |
| `test_evaluate_predictions.R` | 9 | Pearson, Spearman, NDCG@10: perfect/negative correlations, lower-is-better traits, edge cases |
| `test_apply_location_year_scaling.R` | 6 | Z-score standardisation: mean 0 / sd 1 per location-year, NA preservation, zero-variance handling |

### BayesC

| File | Tests | What it covers |
|------|-------|----------------|
| `test_build_fixed_design.R` | 8 | Location + year-within-location design matrix: reference coding, binary output, nested structure |
| `test_build_obs_marker_matrix_and_groups.R` | 7 | Observation-level marker matrix expansion and environment group encoding |
| `test_evaluate_per_location_year.R` | 8 | Index-based per-location-year evaluation: metrics, filtering, NA handling |
| `test_load_and_prepare_data.R` | 12 | Data loading: validation, marker imputation, sample alignment, z-scoring, empty location-year dropping |
| `test_aggregate_helpers.R` | 8 | CSV collection and summary computation for BayesC results aggregation |
| `test_integration_cv_scripts.R` | 18 | End-to-end: CV1, CV2, CV0, CrossLoc scripts produce correct output files with expected structure; fold shuffling between iterations; CV2 stratification; results aggregation |

### RKHS

| File | Tests | What it covers |
|------|-------|----------------|
| `test_rkhs_shared_functions.R` | 14 | Constants match BayesC; design matrix; groups; z-scoring; NDCG; evaluation functions |
| `test_rkhs_build_obs_kernel.R` | 9 | Observation-level kernel expansion: dimensions, symmetry, duplicated genotypes, diagonal values |
| `test_rkhs_kernel_computation.R` | 13 | Gaussian kernel construction: distance matrix properties, kernel symmetry/PSD, bandwidth structure, monomorphic marker removal |
| `test_rkhs_load_and_prepare_data.R` | 11 | Data loading with kernel checkpoint: creation, reproducibility, sample alignment; intentional retention of empty location-years |
| `test_rkhs_fit_and_predict.R` | 8 | BGLR RKHS fitting: groups incompatibility check, prediction columns, variance components, temp file cleanup |
| `test_rkhs_aggregate_helpers.R` | 9 | CSV reading (`read_if_exists`), summary computation, CrossLoc scheme name handling |
| `test_rkhs_integration_cv_scripts.R` | 21 | End-to-end: CV1, CV2, CV0, CrossLoc scripts; seed schemes; fold shuffling between iterations; CV2 stratification; prediction variance; results aggregation |

### GBLUP

| File | Tests | What it covers |
|------|-------|----------------|
| `test_gblup_shared_functions.R` | 17 | Constants match BayesC/RKHS; z-scoring; NDCG; evaluation; summarisation |
| `test_gblup_genotype_alignment.R` | 10 | Factor-level alignment with G matrix: level ordering, extra genotypes, missing genotypes, type coercion, `rowNames` attribute fallback |
| `test_gblup_evaluate_per_location_year.R` | 14 | Join-based per-location-year evaluation: join correctness on (sample.id, location_year), min_genotypes filtering, NA handling, metric verification |
| `test_gblup_fit_and_predict.R` | 8 | ASReml GBLUP fitting: output columns, type checking, location-year filtering, held-out prediction, insufficient data handling |
| `test_gblup_integration_cv_functions.R` | 28 | End-to-end: run_cv1, run_cv2, run_cv0, run_cross_location, run_all_cv_schemes; fold shuffling between iterations; CV2 stratification; seed schemes; cell-level masking; CrossLoc train/predict direction |

### LightGBM (in `models/LightGBM/tests/`)

| File | Tests | What it covers |
|------|-------|----------------|
| `conftest.py` | (fixtures) | Shared fixtures: 100-genotype real data subset, prepared model inputs (PC and kinship), sample predictions DataFrame |
| `test_utils.py` | 25 | Constants, `fetch_model_params`, `calculate_ndcg`, `evaluate_predictions`, `apply_location_year_scaling`, `evaluate_per_location_year`, `summarise_cv_results` |
| `test_prepare_input.py` | 8 | `one_hot_encode`: column creation, binary values, row preservation, column sums |
| `test_tune.py` | 8 | `convert_numpy_types` (JSON serialisation), `tune_trait` (returns params, scoring, skip on small data) |
| `test_run_lightgbm.py` | 36 | `get_feature_columns` (PC/kinship detection, location inclusion/exclusion); `run_cv1`, `run_cv2`, `run_cv0`, `run_cross_location` (output structure, seed schemes, fold shuffling, observation-level masking); evaluation counts (CV2 matches CV1, all location-years in every fold); `run_all_cv_schemes` integration |

## Key test categories

### CV fold assignment and shuffling

The fold-assignment bug (where groups were not shuffled between iterations) is guarded by tests at two levels:

1. **Standalone logic tests** (`test_cv_fold_assignment.R`) — replicate the exact `set.seed()` + `sample()` code from all R model scripts and verify:
   - Different seeds produce different fold compositions
   - Same seed is reproducible
   - CV2 folds are stratified (all location-years in every fold)
   - CV2 folds are balanced within each location-year
   - Unobserved rows stay unassigned

2. **Integration tests** (per-model `test_*integration*.R`) — run 2 iterations of the actual model pipeline and verify from the predictions output that fold membership differs between iterations.

3. **LightGBM** (`test_run_lightgbm.py`) — `test_folds_differ_across_iterations` and `test_cv2_all_location_years_in_every_fold` verify the same properties for `GroupKFold`/`StratifiedKFold`.

### Evaluation metrics

Pearson, Spearman, and NDCG@10 are tested identically across all four model utils (BayesC, RKHS, GBLUP in R; LightGBM in Python). The `lower_is_better` sign flip for DTF, DTH, and plant height is verified in each.

### Per-location-year evaluation

Two variants exist and are tested separately:

- **Index-based** (BayesC, RKHS): `evaluate_per_location_year(observed_data, pred_values, test_indices, trait)` — BGLR returns predictions aligned with input rows
- **Join-based** (GBLUP): `evaluate_per_location_year(observed_data, pred_values, test_ids, trait)` — ASReml returns the full factorial grid, requiring a join on (sample.id, location_year)

### Cross-model consistency

RKHS and GBLUP shared-function tests (`test_rkhs_shared_functions.R`, `test_gblup_shared_functions.R`) explicitly verify that `VALID_TRAITS`, `LOWER_IS_BETTER_TRAITS`, and metric functions match the BayesC reference implementations.
