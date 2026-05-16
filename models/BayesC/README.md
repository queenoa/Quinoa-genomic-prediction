# BayesC Pipeline

BGLR BayesC implementation for the AUSPAK quinoa multi-environment genomic prediction comparison. See the repository [README](../../README.md) for the project overview, trait/environment definitions, CV-scheme designs, and evaluation metrics shared across all four model pipelines.

## Model structure

Defined in `BayesC_utils.R`. Each BGLR fit includes:

- **Fixed effects**: location + year-within-location, via `build_fixed_design()` — approximates an ASReml `fixed = ~ location`, `random = ~ location:year` design, but all fixed since BGLR lacks random effects
- **BayesC marker effects**: SNP-level regression with a spike-and-slab prior (markers can be shrunk to zero)
- **Residual**: per-environment variances via BGLR's `groups` argument (`build_obs_marker_matrix_and_groups()` builds the observation-level group vector alongside the expanded marker matrix)
- **Preprocessing**: z-score standardisation of the trait within each location-year; mean-imputation of missing marker dosages

MCMC settings: 15,000 iterations, 5,000 burn-in, thinning every 5. Overridable via env vars `BAYESC_NITER`, `BAYESC_BURNIN`, `BAYESC_THIN`.

### CrossLoc variant

For cross-location prediction, marker effects are estimated on the source location only, then applied as GEBVs to predict the target location. Implemented in `BayesC_CrossLoc.R`.

## Files

| File | Purpose |
|------|---------|
| `BayesC_utils.R` | Shared utilities: data loading, z-score scaling, marker matrix construction, BGLR model fitting wrapper, evaluation metrics |
| `BayesC_CV1_single_iter.R` | One iteration of CV1 (new genotypes). Called per trait per iteration. |
| `BayesC_CV2_single_iter.R` | One iteration of CV2 (sparse testing). Called per trait per iteration. |
| `BayesC_CV0.R` | CV0 for one trait (deterministic, incremental CSV saving) |
| `BayesC_CrossLoc.R` | CrossLoc for one trait (deterministic, incremental CSV saving) |
| `BayesC_CV0_CrossLoc.R` | Legacy combined CV0 + CrossLoc script (kept for reference) |
| `launch_BayesC_jobs.sh` | Generates and optionally submits SLURM jobs. Supports iteration ranges and selective mode submission. |
| `aggregate_BayesC_results.R` | Collects per-iteration CSVs into combined results, summary statistics, and combined predictions after jobs complete |

## Usage

```bash
# Generate SLURM scripts only
bash launch_BayesC_jobs.sh

# Generate and submit all traits (all CV schemes)
bash launch_BayesC_jobs.sh submit

# Submit a specific iteration range
bash launch_BayesC_jobs.sh submit 1 2
bash launch_BayesC_jobs.sh submit DTF 3 5

# Submit specific CV schemes
bash launch_BayesC_jobs.sh submit cv0
bash launch_BayesC_jobs.sh submit DTF crossloc

# After all CV jobs finish, aggregate results
Rscript aggregate_BayesC_results.R all
```

## Input files

- **Phenotype file**: `../../data/AUSPAK_phenotypes_GP_input.csv` — columns `sample.id`, `location`, `year`, `location_year`, plus one column per trait (per-location means)
- **Marker file**: `../../data/pruned05_AUSPAK_for_bayesC.raw` — PLINK `.raw` format with IID column and SNP dosages (first 6 columns are PLINK metadata)
- **Test data**: `../../data/AUSPAK_test_subset_1k.raw` — 1000-marker subset (551 genotypes) for running tests

## Output files

Per CV iteration / scheme, per trait:
- `cv_results_{SCHEME}_{TRAIT}_..._BayesC.csv` — accuracy metrics per location-year
- `predictions_{SCHEME}_{TRAIT}_..._BayesC.csv` — one row per held-out test accession × location-year (`sample.id`, `location_year`, `location`, `observed`, `predicted`)

After `aggregate_BayesC_results.R`:
- `cv_results_{TRAIT}_BayesC_all_schemes.csv`
- `cv_summary_{TRAIT}_BayesC_all_schemes.csv`
- `predictions_{TRAIT}_BayesC_all_schemes.csv`

## Functions shared with RKHS/GBLUP

Identical-logic helpers in `BayesC_utils.R` (mirrored in `RKHS_utils.R` and `GBLUP_utils.R`):

- `VALID_TRAITS`, `LOWER_IS_BETTER_TRAITS` — trait constants
- `apply_location_year_scaling()` — z-score by location-year
- `build_fixed_design()` — location + year-within-location design matrix (RKHS only; GBLUP uses ASReml formulas instead)
- `calculate_ndcg()` — NDCG@k metric
- `evaluate_predictions()` — Pearson, Spearman, NDCG@10
- `evaluate_per_location_year()` — index-based per-location-year evaluation (BGLR returns row-aligned predictions); GBLUP uses a join-based variant
- `summarise_cv_results()` — mean/SD/min/max per trait × location
