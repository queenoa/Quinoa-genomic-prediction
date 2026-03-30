# GBLUP Multi-Environment Genomic Prediction Pipeline

## Overview

Genomic Best Linear Unbiased Prediction (GBLUP) pipeline for the AUSPAK **quinoa** dataset, using the [ASReml-R](https://vsni.co.uk/software/asreml-r) package. Evaluates prediction accuracy for 7 traits across multiple environments (location-years) using four cross-validation schemes.

Unlike BayesC and RKHS (which use BGLR), GBLUP uses ASReml-R with a pre-computed genomic relationship matrix (G matrix). ASReml supports proper random effects, so location:year is fitted as a random term (shrunk toward zero) rather than approximated as fixed. The model runs significantly faster than BGLR, so all traits and CV schemes are executed in a single monolithic script on a local machine rather than distributed across SLURM jobs.

## GBLUP method

### Genomic relationship matrix

The G matrix inverse (`Ginv_sparse`) is prepared by `G_matrix_GBLUP.R` using the [ASRgenomics](https://cran.r-project.org/package=ASRgenomics) package:

1. Load a pre-computed VanRaden method 1 kinship matrix (computed externally from the full marker set)
2. **Bend** via `ASRgenomics::G.tuneup(G, bend = TRUE)` — adjusts eigenvalues to ensure the G matrix is positive definite (required for stable inversion)
3. **Invert** via `ASRgenomics::G.inverse(Gb, sparse = TRUE)` — computes the inverse and returns it in ASReml's sparse triplet format (a 3-column matrix of row, col, value with genotype IDs stored in the `rowNames` attribute)

ASReml uses `vm(sample.id, Ginv_sparse)` to define the genomic relationship structure. Genotype factor levels are aligned to the G matrix row order via `align_genotypes_to_gmatrix()`, which reads genotype IDs from `attr(Ginv_sparse, "rowNames")`.

### Model structure

Each GBLUP fit (for CV0, CV1, CV2) includes:

- **Fixed effects**: location (main effect)
- **Random effects**: `vm(sample.id, Ginv_sparse)` (genomic BLUPs) + `location:year` (environment-specific deviations, shrunk toward zero)
- **Residual**: `~ units` (single residual variance)
- **Missing data**: `na.action = na.method(y = "include")` — all factor levels remain in the model even when their phenotype is NA, so ASReml can predict for held-out cells
- **Preprocessing**: z-score standardisation of each trait within each location-year

This is the proper mixed-model equivalent of what BayesC/RKHS approximate with all-fixed location + year-within-location design matrices. In GBLUP, `location:year` is random and therefore shrunk toward zero, which is closer to the original ASReml experimental design.

### CrossLoc model variant

For cross-location prediction, the model is simpler:
- **Fixed effects**: intercept only (`trait ~ 1`)
- **Random effects**: `vm(sample.id, Ginv_sparse)` only
- No location or year effects (target location cannot contribute to estimating them)
- Genomic BLUPs (GEBVs) are extracted via `predict(model, classify = "sample.id")` and correlated with observed phenotypes in each target location-year

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

| Scheme | Function | What it tests |
|--------|----------|---------------|
| **CV1** | `run_cv1()` | Predicting **new genotypes** in known environments. 5-fold CV on genotypes, repeated across 15 iterations (seed = 1000+iter). |
| **CV2** | `run_cv2()` | **Sparse testing** — random observation-level cells masked. 5-fold CV stratified by location-year, 15 iterations (seed = 2000+iter). |
| **CV0** | `run_cv0()` | **Leave-one-location-year-out** — entire environments held out. Deterministic, no random folds. |
| **CrossLoc** | `run_cross_location()` | **Cross-location transfer** — train on one location, predict all others via genomic BLUPs (GEBVs). |

All four schemes match the BayesC/RKHS pipelines in design (fold assignment, seed scheme, stratification).

## Evaluation metrics

Predictions are evaluated **per location-year** using:

- **Pearson correlation** — linear predictive accuracy
- **Spearman rank correlation** — rank-based predictive accuracy
- **NDCG@10** — normalized discounted cumulative gain at top 10, for selection ranking quality

Each CV function returns a `predictions` dataframe with one row per test observation, recording `sample.id`, `location_year`, `location`, `observed`, and `predicted` values. These allow post-hoc recomputation of any accuracy metric.

## Files

| File | Purpose |
|------|---------|
| `G_matrix_GBLUP.R` | Prepares `Ginv_sparse`: loads a pre-computed VanRaden kinship matrix, bends it for positive definiteness, computes the sparse inverse, and saves to `Ginv_sparse_GBLUP.RData`. Run once before the CV pipeline. |
| `GBLUP_utils.R` | Shared utilities: constants, z-score scaling, G matrix alignment, ASReml GBLUP fitting wrapper, evaluation metrics, summarisation |
| `GBLUP.R` | Sources utils. Contains all four CV functions (`run_cv1`, `run_cv2`, `run_cv0`, `run_cross_location`), the `run_all_cv_schemes` wrapper, and a usage example. |

## Usage

```bash
cd models/GBLUP

# Step 1: Prepare Ginv (run once)
Rscript G_matrix_GBLUP.R

# Step 2: Run the CV pipeline (redirect output to log file)
Rscript GBLUP.R > GBLUP_run.out 2>&1
```

# Results and predictions are returned in the list:
#   all_cv$cv1$results       — CV1 accuracy metrics
#   all_cv$cv1$predictions   — CV1 individual predictions
#   all_cv$cv1$summary       — CV1 summary statistics
#   all_cv$all_results       — combined results across all schemes
```

## Input data

### Environments

- **Locations (2):** `AUS` (Kununurra, Australia), `PAK` (Faisalabad, Pakistan)
- **Location-years (6):** `AUS_2017`, `AUS_2018`, `AUS_2019`, `PAK_2019`, `PAK_2020`, `PAK_2021`
- Not all genotypes appear in all location-years; not all traits are observed for every genotype × location-year combination (NAs present)

### Raw data files

- **Phenotype file**: data frame with columns `sample.id`, `location`, `year`, `location_year`, and the trait columns
- **Kinship matrix**: VanRaden method 1 kinship matrix (loaded by `G_matrix_GBLUP.R` from an external RData file)
- **G-inverse matrix**: `Ginv_sparse` — ASRgenomics sparse triplet format (3-column matrix of row, col, value) with genotype IDs in the `rowNames` attribute. Produced by `G_matrix_GBLUP.R` and saved to `Ginv_sparse_GBLUP.RData`

## Functions shared with BayesC/RKHS

The following functions in `GBLUP_utils.R` are identical in logic to their counterparts in `BayesC_utils.R` and `RKHS_utils.R`:

- `VALID_TRAITS`, `LOWER_IS_BETTER_TRAITS` — trait constants
- `apply_location_year_scaling()` — z-score standardisation by location-year
- `calculate_ndcg()` — NDCG@k metric
- `evaluate_predictions()` — Pearson, Spearman, NDCG@10 evaluation
- `summarise_cv_results()` — mean/SD/min/max summary per location

## Why `evaluate_per_location_year()` differs between GBLUP and BGLR pipelines

The BayesC and RKHS pipelines (both using BGLR) use an **index-based** `evaluate_per_location_year()`:

```r
evaluate_per_location_year(observed_data, pred_values, test_indices, trait)
```

BGLR returns predictions as a data frame with one row per observation in the same order as the input data. The function indexes directly into both `observed_data` and `pred_values` by row position (`test_indices`), so evaluation is a simple positional lookup.

GBLUP uses a **join-based** `evaluate_per_location_year()`:

```r
evaluate_per_location_year(observed_data, pred_values, test_ids, trait)
```

ASReml's `predict(model, classify = "sample.id:location:year")` returns the full factorial grid of all genotype x location x year combinations (including ones that don't exist in the data). The output is not row-aligned with the input data, so evaluation requires joining on `(sample.id, location_year)` to match predictions to observations.

Both approaches produce the same result — Pearson, Spearman, and NDCG@10 per location-year for the held-out observations. The difference is purely mechanical, driven by how each package returns predictions.

**Important for CV2**: The join-based function uses genotype IDs (`test_ids`), not row indices. In CV2, masking is at the cell level — the same genotype may be masked in one location-year but observed in another. Using `evaluate_per_location_year()` with genotype IDs would match unmasked observations too, causing data leakage. For this reason, CV2 uses inline cell-level evaluation (tracking masked cells explicitly and joining by `(sample.id, location_year)`) rather than calling `evaluate_per_location_year()`.
