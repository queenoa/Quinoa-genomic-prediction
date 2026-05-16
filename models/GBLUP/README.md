# GBLUP Pipeline

ASReml-R GBLUP implementation for the AUSPAK quinoa multi-environment genomic prediction comparison. See the repository [README](../../README.md) for the project overview, trait/environment definitions, CV-scheme designs, and evaluation metrics shared across all four pipelines.

Two model variants are run:

- **Global model** — fitted on data pooled across both locations. Primary pipeline (`GBLUP.R` + `GBLUP_utils.R`). Used for all four CV schemes (CV1, CV2, CV0, CrossLoc).
- **Per-location model** — fitted independently on each location's subset (`GBLUP_per_location_CV1.R`). Used to benchmark within-location vs. global prediction accuracy and to estimate location-specific genetic variance components.

GBLUP uses ASReml-R with a pre-computed genomic relationship matrix. Runs locally (no SLURM).

## Genomic relationship matrix

The G matrix inverse (`Ginv_sparse`) is prepared by `G_matrix_GBLUP.R` using [ASRgenomics](https://cran.r-project.org/package=ASRgenomics):

1. Load a pre-computed VanRaden method 1 kinship matrix
2. **Bend** via `ASRgenomics::G.tuneup(G, bend = TRUE)` — adjusts eigenvalues to ensure positive definiteness
3. **Invert** via `ASRgenomics::G.inverse(Gb, sparse = TRUE)` — sparse triplet format (row, col, value) with genotype IDs in the `rowNames` attribute

ASReml uses `vm(sample.id, Ginv_sparse)` to define the genomic relationship. Genotype factor levels are aligned to G matrix row order via `align_genotypes_to_gmatrix()`, which reads genotype IDs from `attr(Ginv_sparse, "rowNames")`.

## Global model structure

Pooled across both locations. Used for CV1 and CV2. Implemented in `fit_gblup_and_predict()` (`GBLUP_utils.R`):

- **Fixed**: `trait ~ location`
- **Random**: `vm(sample.id, Ginv_sparse) + at(location):year` — genomic BLUPs plus **location-specific year effects with separate variance components per location**
- **Residual**: `dsum(~ units | location)` — **heterogeneous residual variance per location**
- **Missing data**: `na.action = na.method(y = "include")` — all factor levels remain even when phenotype is NA, so ASReml predicts held-out cells
- **Preprocessing**: z-score standardisation per location-year
- **Data ordering**: rows are sorted by `location` inside `fit_gblup_and_predict()` (required by `dsum()`)

`fit_gblup_and_predict()` returns predictions and variance components (`summary(model)$varcomp`). The CV runners annotate the varcomps with per-fit metadata (`trait`, `iteration`, `fold`, `seed`, `cv_scheme`) and propagate them to `varcomps_*` CSVs.

### CV0 variant

CV0 uses a separate, simpler asreml fit defined inline in `run_cv0()` (`GBLUP.R`):

- **Fixed**: `trait ~ location`
- **Random**: `vm(sample.id, Ginv_sparse) + location:year` — homogeneous year-within-location variance
- **Residual**: `~ units` — **homogeneous residual variance**

The heterogeneous spec (`dsum(~ units | location)` + `at(location):year`) fails for the AUS traits that have only two of three years available when the held-out location-year is removed. The homogeneous structure pools across locations and years and avoids these failures.

### CrossLoc variant

For cross-location prediction (`run_cross_location`):

- **Fixed**: intercept only (`trait ~ 1`)
- **Random**: `vm(sample.id, Ginv_sparse)` only
- **Residual**: `~ units` (within a single location, `dsum(~ units | location)` and `at(location):year` collapse to homogeneous)
- GEBVs are extracted via `predict(model, classify = "sample.id")` and correlated with observed phenotypes in each target location-year. Variance components from each within-location fit are captured and saved.

## Per-location models

Implemented in `GBLUP_per_location_CV1.R`. Fitted independently on each location's data subset. Two analyses use the same model:

1. **Variance-component estimation** (`fit_single_location_gblup` / `run_per_location_variances`) — fit once per `(trait, location)` on the full data, extract `genetic_variance`, `residual_variance`, and `h2 = Vg / (Vg + Ve)` from `summary(model)$varcomp`. Saved to `per_location_variances_GBLUP.csv`.
2. **CV1 within location** (`fit_gblup_per_location` / `run_cv1_per_location`) — 5-fold CV on genotypes within each location, 15 iterations, seed = `1000 + iter`. Matches the global CV1 fold-assignment scheme. Produces per `(trait, location)` and combined-per-trait CSVs.

Per-location model structure (both analyses):

- **Fixed**: intercept only (`trait ~ 1`) — single location, no location effect to fit
- **Random**: `vm(sample.id, Ginv_sparse) + year` — homogeneous year variance within the location; falls back to `vm(sample.id, Ginv_sparse)` only when the location has a single year
- **Residual**: `~ units` (homogeneous)
- **Preprocessing**: z-score per location-year (same as global)

**Why homogeneous residuals at the per-location scale?** Heterogeneous residual variances across years (`dsum(~ units | year)`) were not modelled in the per-location models due to frequent convergence failures, given the limited number of years per location (2–3), the absence of within-year replication, and the presence of `year` as a random effect.

## CV-scheme → function/script map

| Scheme | Function | Script |
|--------|----------|--------|
| CV1 | `run_cv1()` | `GBLUP.R` |
| CV1 per-location | `run_cv1_per_location()` | `GBLUP_per_location_CV1.R` |
| CV2 | `run_cv2()` | `GBLUP.R` |
| CV0 | `run_cv0()` | `GBLUP.R` |
| CrossLoc | `run_cross_location()` | `GBLUP.R` |

Fold assignment, seed scheme, and stratification match the BayesC/RKHS pipelines.

## Files

| File | Purpose |
|------|---------|
| `G_matrix_GBLUP.R` | Prepares `Ginv_sparse_GBLUP.RData` from a pre-computed VanRaden kinship matrix (bend + sparse invert). Run once. |
| `GBLUP_utils.R` | Shared utilities for the global pipeline: constants, z-score scaling, G-alignment, `fit_gblup_and_predict()` (heterogeneous global model used by CV1/CV2, returning predictions + varcomps), join-based per-location-year evaluation, summarisation. |
| `GBLUP.R` | Global pipeline. All four CV functions (`run_cv1`, `run_cv2`, `run_cv0`, `run_cross_location`), the `run_all_cv_schemes` wrapper, and CSV output blocks (results, predictions, varcomps). `run_cv0` and `run_cross_location` define their asreml fits inline (homogeneous residuals); `run_cv1` and `run_cv2` go through `fit_gblup_and_predict()`. |
| `GBLUP_per_location_CV1.R` | Per-location pipeline. Fits the within-location GBLUP per (trait, location): (a) once on full data to extract variance components → `per_location_variances_GBLUP.csv`; (b) under a 5-fold × 15-iter CV1 design → per-location CV1 CSVs. |

## Usage

```bash
cd models/GBLUP

# Step 1: Prepare Ginv (run once)
Rscript G_matrix_GBLUP.R

# Step 2: Global CV pipeline (all four schemes)
Rscript GBLUP.R > GBLUP_run.out 2>&1

# Step 3: Per-location CV1 + per-location variance components
Rscript GBLUP_per_location_CV1.R > GBLUP_per_location_run.out 2>&1
```

In-memory return from `run_all_cv_schemes(...)` (the global pipeline):

- `$cv1$results / $predictions / $varcomps / $summary` — and likewise `$cv2`, `$cv0`, `$cross_loc`
- `$all_results` — combined results across all four schemes

## Input files

- **Phenotype file**: `../../data/AUSPAK_phenotypes_GP_input.csv` — columns `sample.id`, `location`, `year`, `location_year`, plus one column per trait (per-location means)
- **Kinship matrix**: VanRaden method 1 (loaded by `G_matrix_GBLUP.R` from an external RData file)
- **G-inverse matrix**: `Ginv_sparse_GBLUP.RData` — ASRgenomics sparse triplet format with genotype IDs in the `rowNames` attribute; produced by `G_matrix_GBLUP.R`

## Output files

Per scheme, per trait (from `GBLUP.R`):
- `cv_results_{SCHEME}_{TRAIT}_GBLUP.csv` — accuracy metrics per location-year
- `predictions_{SCHEME}_{TRAIT}_GBLUP.csv` — individual held-out predictions
- `varcomps_{SCHEME}_{TRAIT}_GBLUP.csv` — variance components (one row per `summary(model)$varcomp` entry, annotated with `trait`, `iteration`, `fold`, `seed`, `cv_scheme`)

Combined across schemes, per trait (from `GBLUP.R`):
- `cv_results_{TRAIT}_GBLUP_all_schemes.csv`
- `cv_summary_{TRAIT}_GBLUP_all_schemes.csv` — mean / SD / min / max per scheme × location
- `predictions_{TRAIT}_GBLUP_all_schemes.csv`
- `varcomps_{TRAIT}_GBLUP_all_schemes.csv`

Per-location pipeline (from `GBLUP_per_location_CV1.R`):
- `per_location_variances_GBLUP.csv` — full-data variance components per `(trait, location)`: `n_obs`, `genetic_variance`, `residual_variance`, `h2`
- `cv_results_CV1_{LOC}_{TRAIT}_GBLUP.csv` and `predictions_CV1_{LOC}_{TRAIT}_GBLUP.csv` — per-location CV1 outputs
- `cv_results_{TRAIT}_GBLUP_CV1_per_location.csv` and `predictions_{TRAIT}_GBLUP_CV1_per_location.csv` — combined across locations per trait
- `cv_summary_GBLUP_CV1_per_location.csv` — summary across all traits and locations

## Functions shared with BayesC/RKHS

Identical-logic helpers in `GBLUP_utils.R` (mirrored in `BayesC_utils.R` and `RKHS_utils.R`):

- `VALID_TRAITS`, `LOWER_IS_BETTER_TRAITS`
- `apply_location_year_scaling()` — z-score by location-year
- `calculate_ndcg()` — NDCG@k
- `evaluate_predictions()` — Pearson, Spearman, NDCG@10
- `summarise_cv_results()` — mean / SD / min / max per trait × location

## Why `evaluate_per_location_year()` differs between GBLUP and BGLR pipelines

BayesC and RKHS (both BGLR) use an **index-based** `evaluate_per_location_year(observed_data, pred_values, test_indices, trait)`. BGLR returns predictions row-aligned with the input data, so evaluation is a positional lookup.

GBLUP uses a **join-based** `evaluate_per_location_year(observed_data, pred_values, test_ids, trait)`. ASReml's `predict(model, classify = "sample.id:location:year")` returns the full factorial grid of all genotype × location × year combinations (including non-existent ones), so evaluation requires joining on `(sample.id, location_year)` to match predictions back to observations.

**CV2 caveat**: the join-based function uses genotype IDs (`test_ids`), not row indices. In CV2, masking is at the cell level — the same genotype may be masked in one location-year but observed in another. Using the join-based helper directly would match unmasked observations too, causing data leakage. CV2 therefore tracks masked cells explicitly and joins by `(sample.id, location_year)` inline rather than calling `evaluate_per_location_year()`.
