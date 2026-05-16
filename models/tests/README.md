# Test Suite for AUSPAK Genomic Selection Models

## Overview

Shared test directory covering all four model pipelines: **BayesC**, **RKHS**, **GBLUP** (R, testthat), and **LightGBM** (Python, pytest). Tests range from lightweight unit tests of individual utility functions to end-to-end integration tests that fit models on a 1000-marker test subset.

All test files live together in `models/tests/`. Python (LightGBM) and R (BayesC/RKHS/GBLUP) tests coexist; pytest and testthat each pick up their respective files.

## Running tests

### R tests (BayesC, RKHS, GBLUP)

```bash
cd models/tests

# Run all R tests
Rscript -e 'testthat::test_dir(".")'

# Run a single test file
Rscript -e 'testthat::test_file("test_cv_fold_assignment.R")'

# Run only tests matching a pattern
Rscript -e 'testthat::test_dir(".", filter = "gblup")'
```

**Dependencies:** `testthat`, `dplyr`, `data.table`, `BGLR` (BayesC/RKHS integration tests), `asreml` + `ASRgenomics` (GBLUP integration tests).

**Test data:** `../../data/AUSPAK_test_subset_1k.raw` (1000-marker subset, 551 genotypes) and `../../data/AUSPAK_phenotypes_GP_input.csv` (per-location-year means — not BLUEs; the largely unreplicated trial structure made reliable within-environment BLUEs infeasible).

### Python tests (LightGBM)

```bash
cd models/tests
conda activate ../LightGBM/env-ML
pytest -v
```

**Dependencies:** `pytest`, `scikit-learn`, `lightgbm`, `pandas`, `numpy`.

The LightGBM `conftest.py` adds `../LightGBM/` to `sys.path` so test files can `from run_LightGBM import ...`.

## Test fixtures

`fixtures/` holds 100-genotype subsets of the real AUSPAK data, used by the Python tests:

- `pheno_subset.csv` — phenotype subset with the production trait columns (`DTF`, `DTH`, `PtHt`, `PcleLng`, `SdLen`, `TGW`, `SdW_z` — per-location means)
- `kinship_subset.csv` — 100×100 kinship matrix subset

Regenerate via `python models/tests/create_test_fixtures.py` (needs the env-ML conda environment).

## Test speed tiers

| Tier | Files | Time | What runs |
|------|-------|------|-----------|
| **Instant** (<1s) | `test_cv_fold_assignment.R`, `test_calculate_ndcg.R`, `test_evaluate_predictions.R`, `test_apply_location_year_scaling.R`, `test_build_fixed_design.R`, `test_build_obs_marker_matrix_and_groups.R` | <1s each | Pure R logic, no model fitting |
| **Fast** (1-10s) | `test_evaluate_per_location_year.R`, `test_load_and_prepare_data.R`, `test_aggregate_helpers.R`, `test_gblup_shared_functions.R`, `test_gblup_genotype_alignment.R`, `test_gblup_evaluate_per_location_year.R`, `test_rkhs_shared_functions.R`, `test_rkhs_build_obs_kernel.R`, `test_rkhs_aggregate_helpers.R`, `test_rkhs_kernel_computation.R` | 1-10s each | Data loading, matrix ops, no model fitting |
| **Slow** (1-10 min) | `test_integration_cv_scripts.R`, `test_rkhs_integration_cv_scripts.R`, `test_gblup_integration_cv_functions.R`, `test_gblup_per_location.R`, `test_rkhs_fit_and_predict.R`, `test_gblup_fit_and_predict.R`, `test_rkhs_load_and_prepare_data.R` | 1-10 min each | BGLR/ASReml model fitting with minimal MCMC |

To run only the fast tests (no model fitting):

```bash
Rscript -e 'testthat::test_dir(".", filter = "^(?!.*integration|.*fit_and_predict|.*per_location|.*rkhs_load)")'
```

## Trait names

All tests use the production trait names without suffix: `DTF`, `DTH`, `PtHt`, `PcleLng`, `SdLen`, `TGW`, `SdW_z`. The lower-is-better set is `DTF`, `DTH`, `PtHt`.

## Integration tests — directory layout

The BayesC and RKHS `*_utils.R` files load the phenotype CSV from the default relative path `"../AUSPAK_phenotypes_GP_input.csv"`. The integration tests therefore set up a parent-and-work layout in `tempdir()` so the default path resolves correctly:

```
PARENT_DIR/                              (== OUTDIR/..)
├── AUSPAK_phenotypes_GP_input.csv       (matches default)
└── work/                                (== OUTDIR; CWD when scripts run)
    └── (utils file, output CSVs)
```

The test runs the CV script with CWD set to `OUTDIR`, so `load_and_prepare_data()`'s default `"../AUSPAK_phenotypes_GP_input.csv"` finds the copied phenotype CSV in `PARENT_DIR`.

## Test files by model

### Cross-model (shared logic)

| File | What it covers |
|------|----------------|
| `test_cv_fold_assignment.R` | CV1 and CV2 fold creation logic shared by all R models: shuffling across iterations, reproducibility, balance, stratification, unobserved row handling |
| `test_calculate_ndcg.R` | NDCG@k: perfect/worst rankings, lower-is-better, ties, negatives, edge cases |
| `test_evaluate_predictions.R` | Pearson, Spearman, NDCG@10: perfect/negative correlations, lower-is-better traits, edge cases |
| `test_apply_location_year_scaling.R` | Z-score standardisation: mean 0 / sd 1 per location-year, NA preservation, zero-variance handling |

### BayesC

| File | What it covers |
|------|----------------|
| `test_build_fixed_design.R` | Location + year-within-location design matrix: reference coding, binary output, nested structure |
| `test_build_obs_marker_matrix_and_groups.R` | Observation-level marker matrix expansion and environment group encoding |
| `test_evaluate_per_location_year.R` | Index-based per-location-year evaluation: metrics, filtering, NA handling |
| `test_load_and_prepare_data.R` | Data loading: validation, marker imputation, sample alignment, marker centring, z-scoring, empty location-year dropping |
| `test_aggregate_helpers.R` | CSV collection and summary computation for BayesC results aggregation |
| `test_integration_cv_scripts.R` | End-to-end: CV1, CV2, **CV0** (`BayesC_CV0.R`), and **CrossLoc** (`BayesC_CrossLoc.R`) scripts produce correct output files with expected structure; fold shuffling between iterations; CV2 stratification; results aggregation |

### RKHS

| File | What it covers |
|------|----------------|
| `test_rkhs_shared_functions.R` | Constants match BayesC; design matrix; groups; z-scoring; NDCG; evaluation functions |
| `test_rkhs_build_obs_kernel.R` | Observation-level kernel expansion: dimensions, symmetry, duplicated genotypes, diagonal values |
| `test_rkhs_kernel_computation.R` | Gaussian kernel construction: distance matrix properties, kernel symmetry/PSD, bandwidth structure, monomorphic marker removal |
| `test_rkhs_load_and_prepare_data.R` | Data loading with kernel checkpoint: creation, reproducibility, sample alignment; intentional retention of empty location-years |
| `test_rkhs_fit_and_predict.R` | BGLR RKHS fitting: groups incompatibility check, prediction columns, variance components, temp file cleanup |
| `test_rkhs_aggregate_helpers.R` | CSV reading (`read_if_exists`), summary computation, CrossLoc scheme name handling |
| `test_rkhs_integration_cv_scripts.R` | End-to-end: CV1, CV2, CV0, CrossLoc scripts; seed schemes; fold shuffling between iterations; CV2 stratification; prediction variance; results aggregation |

### GBLUP

| File | What it covers |
|------|----------------|
| `test_gblup_shared_functions.R` | Constants match BayesC/RKHS; z-scoring; NDCG; evaluation; summarisation |
| `test_gblup_genotype_alignment.R` | Factor-level alignment with G matrix: level ordering, extra genotypes, missing genotypes, type coercion, `rowNames` attribute fallback |
| `test_gblup_evaluate_per_location_year.R` | Join-based per-location-year evaluation: join correctness on (sample.id, location_year), min_genotypes filtering, NA handling, metric verification |
| `test_gblup_fit_and_predict.R` | ASReml GBLUP fitting (heterogeneous global model used by CV1/CV2): `pred_values` / `varcomp` list shape; held-out predictions; insufficient-data handling; varcomp annotation; resilience to unsorted input rows (required by `dsum`) |
| `test_gblup_integration_cv_functions.R` | End-to-end: `run_cv1`, `run_cv2`, `run_cv0`, `run_cross_location`, `run_all_cv_schemes`; fold shuffling between iterations; CV2 stratification; seed schemes; cell-level masking; CrossLoc train/predict direction; **per-scheme `varcomps` capture with correct annotation (`iteration`/`fold`/`seed`/`cv_scheme` for CV1/CV2, `held_out_location_year` for CV0, `train_location` for CrossLoc); CV0 uses the inline homogeneous fit (`location:year` + `~ units`), not the heterogeneous `dsum`/`at(location):year` structure** |
| `test_gblup_per_location.R` | **`GBLUP_per_location_CV1.R` functions**: `fit_gblup_per_location` (single-location CV fit, predictions stay within focal location), `fit_single_location_gblup` (full-data variance components — Vg, Ve, h2), `run_per_location_variances` (one row per trait × location), `run_cv1_per_location` (CV1 design, train_location tagging, CV1 seed scheme) |

`GBLUP.R` and `GBLUP_per_location_CV1.R` carry a `sys.nframe() == 0` guard so the run-pipeline block at the end of each file is skipped when sourced from inside tests; only the function definitions are loaded.

### LightGBM

| File | What it covers |
|------|----------------|
| `conftest.py` | Shared fixtures: 100-genotype real data subset (`pheno_subset`, `kinship_subset`), prepared kinship model input (`model_input_kinship`), sample predictions DataFrame (`predictions_df`). Adds `../LightGBM/` to `sys.path` |
| `test_utils.py` | Constants (`VALID_TRAITS`, `LOWER_IS_BETTER_TRAITS`, `DEFAULT_PARAMS`), `fetch_model_params`, `calculate_ndcg`, `evaluate_predictions`, `apply_location_year_scaling`, `evaluate_per_location_year`, `summarise_cv_results` |
| `test_prepare_input.py` | `one_hot_encode`: column creation, binary values, row preservation, column sums (kinship + pheno merge mirroring the production pipeline) |
| `test_tune.py` | `convert_numpy_types` (JSON serialisation), `tune_trait` (returns params, scoring, skip on small data) |
| `test_run_lightgbm.py` | `get_feature_columns` (kinship detection, location inclusion/exclusion, **PC columns explicitly NOT detected** since they were dropped from the pipeline in commit c97e97a); `run_cv1`, `run_cv2`, `run_cv0`, `run_cross_location` (output structure, seed schemes, fold shuffling, observation-level masking); evaluation counts (CV2 matches CV1, all location-years in every fold); `run_all_cv_schemes` integration |

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

### GBLUP variance components

The heterogeneous-residual global model (CV1/CV2 via `fit_gblup_and_predict()`) returns variance components alongside predictions. The CV runners annotate each varcomp row with per-fit metadata so the `varcomps_*.csv` outputs can be analysed across folds/iterations. Verified by `test_gblup_fit_and_predict.R` (single-fit return shape) and `test_gblup_integration_cv_functions.R` (per-scheme annotation correctness). CV0 fits its own asreml model inline with **homogeneous residuals** (`residual = ~ units`, `random = ~ vm(sample.id, Ginv_sparse) + location:year`) because the heterogeneous `dsum`/`at(location):year` structure fails for AUS traits with only two of three years when a whole location-year is held out; `test_gblup_integration_cv_functions.R` asserts CV0 varcomps contain no per-location `AUS!…`/`PAK!…` residuals nor `at(location):year` random components.
