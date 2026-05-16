# LightGBM Pipeline

Gradient-boosting (LightGBM) implementation for the AUSPAK quinoa multi-environment genomic prediction comparison. See the repository [README](../../README.md) for the project overview, trait/environment definitions, CV-scheme designs, and evaluation metrics shared across all four model pipelines.

Runs locally (no SLURM). All traits and CV schemes are executed in a single script.

## Features

Genetic features:
- **Kinship matrix**: VanRaden kinship loaded from `kinship_matrix_VanRaden_auspak_maxmissing20.csv`. Columns are prefixed `K_` (e.g. `K_S3H3_batch1`). Each row's kinship features represent its genetic relatedness to all other genotypes.

Environment features:
- One-hot-encoded `location` and `location_year`. These are the LightGBM analogue of the fixed effects in the BayesC/RKHS/GBLUP models.

`get_feature_columns()` in `run_LightGBM.py` selects all `K_`-prefixed columns plus the environment encodings. The CrossLoc scheme uses kinship features only (no location/location_year encoding) — see below.

## Model parameters

Default (pre-tuning) parameters:

```python
DEFAULT_PARAMS = {'max_depth': 3, 'learning_rate': 0.05, 'n_estimators': 500}
```

Per-trait tuned parameters are produced by `tune_LightGBM.py` and written to `tuned_params.json`. `run_LightGBM.py` loads them automatically when present; otherwise it falls back to `DEFAULT_PARAMS`.

## Hyperparameter tuning

`tune_LightGBM.py` runs a single `RandomizedSearchCV` per trait on the full dataset with 5-fold `GroupKFold` (grouped by genotype) to select reasonable hyperparameters. 50 random combinations are explored over:

| Parameter | Search values |
|-----------|---------------|
| `n_estimators` | 100, 200, 500, 1000 |
| `max_depth` | 2, 3, 4, 5, 6, -1 |
| `learning_rate` | 0.01, 0.03, 0.05, 0.1, 0.2 |
| `num_leaves` | 15, 31, 63, 127 |
| `min_child_samples` | 5, 10, 20, 50 |
| `subsample` | 0.6, 0.7, 0.8, 0.9, 1.0 |
| `colsample_bytree` | 0.6, 0.7, 0.8, 0.9, 1.0 |
| `reg_alpha` | 0, 0.01, 0.1, 1.0 |
| `reg_lambda` | 0, 0.01, 0.1, 1.0 |

Scoring: `neg_mean_squared_error`. Best parameters per trait are saved to `tuned_params.json`. One-shot step — run once before the CV pipeline.

## Preprocessing

Z-score standardisation of each trait within each location-year, applied once upfront before all CV schemes. Identical logic to the BayesC/RKHS/GBLUP `apply_location_year_scaling()`.

## CrossLoc feature handling

For cross-location prediction, only kinship features are used (no location or location_year encoding). This matches GBLUP's `trait ~ 1 + G` cross-location model — prediction comes purely from genetic relationships captured in the kinship matrix.

## Per-location-year evaluation

Minimum 10 genotypes per location-year to compute metrics. When multiple observations per genotype exist within a location-year, predictions are averaged by genotype before evaluation.

## Files

| File | Purpose |
|------|---------|
| `Prepare_input_data.py` | Data preparation: loads phenotypes and kinship matrix; produces `model_inputs/model_input.pkl` with one-hot-encoded `location` and `location_year`. Run once. |
| `tune_LightGBM.py` | Hyperparameter tuning: runs `RandomizedSearchCV` per trait (50 iterations, 5-fold `GroupKFold`), saves best params to `tuned_params.json`. Run once. |
| `LightGBM_utils.py` | Shared utilities: constants, `DEFAULT_PARAMS`, `fetch_model_params()`, z-score scaling, evaluation metrics (Pearson, Spearman, NDCG@10), per-location-year evaluation, summarisation |
| `run_LightGBM.py` | Main pipeline: loads tuned params (or defaults), contains all four CV functions (`run_cv1`, `run_cv2`, `run_cv0`, `run_cross_location`), the `run_all_cv_schemes` wrapper, and CSV output |
| `environment_ml.yml` | Conda environment specification |
| `README_env_ml.txt` | Environment setup notes |

## Usage

```bash
cd models/LightGBM

# Step 1: Conda environment (once)
conda env create --file environment_ml.yml --prefix ./env-ML
conda activate ./env-ML

# Step 2: Prepare input data (once → model_inputs/model_input.pkl)
python Prepare_input_data.py

# Step 3: Tune hyperparameters (once → tuned_params.json)
python tune_LightGBM.py 2>&1 | tee tune_LightGBM.out

# Step 4: Run the CV pipeline
python run_LightGBM.py 2>&1 | tee run_LightGBM.out
```

## Input files

- **Phenotype file**: `../../data/AUSPAK_phenotypes_GP_input.csv` — columns `location_year`, `location`, `year`, `sample.id`, plus one column per trait (`DTF`, `DTH`, `PtHt`, `PcleLng`, `SdLen`, `TGW`, `SdW_z`) containing per-location means
- **Kinship file**: `../../data/kinship_matrix_VanRaden_auspak_maxmissing20.csv` — VanRaden kinship matrix with sample IDs as row index and column names

### Prepared pickle (produced by `Prepare_input_data.py`)

`model_inputs/model_input.pkl` has one row per genotype × location-year observation, with column groups:

- **Genetic features**: `K_*` (kinship)
- **Trait columns**: `DTF`, `DTH`, `PtHt`, `PcleLng`, `SdLen`, `TGW`, `SdW_z` (per-location means; may contain NAs)
- **Metadata**: `sample.id`, `location_year`, `location` (retained for CV splitting/evaluation; never used as features)
- **One-hot environment features**: `location_*` and `location_year_*` (produced by `one_hot_encode()`), used as features by `get_feature_columns()`

### Tuned parameters

`tuned_params.json` — per-trait LightGBM hyperparameters from `tune_LightGBM.py`.

## Output files

Per scheme, per trait:
- `cv_results_{SCHEME}_{TRAIT}_LightGBM.csv` — accuracy metrics per location-year
- `predictions_{SCHEME}_{TRAIT}_LightGBM.csv` — individual predictions

Combined all-schemes per trait:
- `cv_results_{TRAIT}_LightGBM_all_schemes.csv` — all CV results
- `cv_summary_{TRAIT}_LightGBM_all_schemes.csv` — summary statistics (mean/SD/min/max)
- `predictions_{TRAIT}_LightGBM_all_schemes.csv` — all predictions

## Functions shared with BayesC/RKHS/GBLUP

The following in `LightGBM_utils.py` are identical in logic to their R counterparts:

- `VALID_TRAITS`, `LOWER_IS_BETTER_TRAITS`
- `apply_location_year_scaling()` — z-score by location-year
- `calculate_ndcg()` — NDCG@k
- `evaluate_predictions()` — Pearson, Spearman, NDCG@10
- `summarise_cv_results()` — mean/SD/min/max per trait × location

## LightGBM-specific utilities

In `LightGBM_utils.py`:

- `DEFAULT_PARAMS` — pre-tuning fallback hyperparameters (`max_depth: 3`, `learning_rate: 0.05`, `n_estimators: 500`)
- `fetch_model_params()` — looks up per-trait hyperparameters from a tuned-params dict; falls back to `DEFAULT_PARAMS` if the dict is flat or `None`
