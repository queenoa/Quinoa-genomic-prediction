# LightGBM Multi-Environment Genomic Prediction Pipeline

## Overview

Gradient boosting (LightGBM) pipeline for the AUSPAK quinoa dataset. Evaluates prediction accuracy for 7 traits across multiple environments (location-years) using four cross-validation schemes.

Unlike BayesC, RKHS (which use BGLR), and GBLUP (which uses ASReml), LightGBM supports three alternative genetic feature sets: all PCs, 25 PCs, or the kinship matrix. Kinship matrix was chosen for the final model. Environment effects are captured through one-hot encoded location and location-year variables. The model runs fast enough that all traits and CV schemes are executed in a single script (like GBLUP), with no SLURM job distribution needed.

## LightGBM method

### Feature construction

Three alternative genetic feature sets are supported, each produced as a separate pickle by `Prepare_input_data.py`:

1. **All PCs (551)**: Principal components computed from the full marker set, loaded from `AUSPAK_PCs_all.csv`. Columns prefixed `PC`. This is the default.
2. **25 PCs**: First 25 principal components (PC1–PC25) — a reduced-dimensionality alternative.
3. **Kinship matrix**: VanRaden kinship matrix loaded from `kinship_matrix_VanRaden_auspak_maxmissing20.csv`. Columns are prefixed `K_` (e.g. `K_S3H3_batch1`) to distinguish them from PCs. Each row's kinship features represent its genetic relatedness to all other genotypes.

Environment features (shared across all three):
- One-hot encoded `location` and `location_year` columns. Location captures main location effects; location_year captures year-within-location effects. These are analogous to the fixed effects in the BayesC/RKHS/GBLUP models.

`get_feature_columns()` in `run_LightGBM.py` auto-detects genetic features by prefix (`PC` or `K_`), so switching feature set only requires changing which pickle is loaded.

### Model parameters

Default (pre-tuning) parameters:

```python
DEFAULT_PARAMS = {'max_depth': 3, 'learning_rate': 0.05, 'n_estimators': 500}
```

Per-trait tuned parameters are obtained via `tune_LightGBM.py` (see **Hyperparameter tuning** below) and saved to `tuned_params.json`. When present, `run_LightGBM.py` loads them automatically; otherwise falls back to defaults.

### Hyperparameter tuning

`tune_LightGBM.py` runs a single `RandomizedSearchCV` per trait on the full dataset with 5-fold `GroupKFold` (grouped by genotype) to select reasonable hyperparameters. The search explores 50 random combinations over:

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

Scoring: `neg_mean_squared_error`. Best parameters per trait are saved to `tuned_params.json`. This is a one-shot step — run once before the CV pipeline.

### Preprocessing

Z-score standardisation of each trait within each location-year, applied once upfront before all CV schemes. Identical logic to BayesC/RKHS/GBLUP `apply_location_year_scaling()`.

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

Traits where lower is better (DTF, DTH, plant height) have NDCG@10 sign-flipped accordingly.

## Cross-validation schemes

| Scheme | Function | What it tests | Features used |
|--------|----------|---------------|---------------|
| **CV1** | `run_cv1()` | Predicting **new genotypes** in known environments. 5-fold CV on genotypes, 15 iterations (seed = 1000+iter). | genetic features + location + location_year |
| **CV2** | `run_cv2()` | **Sparse testing** — random observation-level cells masked. 5-fold CV stratified by location-year, 15 iterations (seed = 2000+iter). | genetic features + location + location_year |
| **CV0** | `run_cv0()` | **Leave-one-location-year-out** — entire environments held out. Deterministic, no random folds. | genetic features + location + location_year |
| **CrossLoc** | `run_cross_location()` | **Cross-location transfer** — train on one location, predict all others. | genetic features only (no location/location_year encoding) |

All four schemes match the BayesC/RKHS/GBLUP pipelines in design (fold assignment, seed scheme, stratification).

### CrossLoc feature handling

For cross-location prediction, only genetic features are used (no location or location_year encoding). This matches GBLUP's approach of using only `trait ~ 1 + G` with no location effects — prediction comes purely from genetic relationships captured in the PCs or kinship matrix.

## Evaluation metrics

Predictions are evaluated **per location-year** using:

- **Pearson correlation** — linear predictive accuracy
- **Spearman rank correlation** — rank-based predictive accuracy
- **NDCG@10** — normalised discounted cumulative gain at top 10, for selection ranking quality

Minimum 10 genotypes per location-year to compute metrics. When multiple observations per genotype exist within a location-year, predictions are averaged by genotype before evaluation.

## Files

| File | Purpose |
|------|---------|
| `Prepare_input_data.py` | Data preparation: loads phenotypes, PCA data, and kinship matrix; produces three pickles (all PCs, 25 PCs, kinship) with one-hot encoded location and location_year. Run once before the CV pipeline. |
| `tune_LightGBM.py` | Hyperparameter tuning: runs RandomizedSearchCV per trait (50 iterations, 5-fold GroupKFold), saves best params to `tuned_params.json`. Run once before the CV pipeline. |
| `LightGBM_utils.py` | Shared utilities: constants, default params, param resolver, z-score scaling, evaluation metrics (Pearson, Spearman, NDCG@10), per-location-year evaluation, summarisation |
| `run_LightGBM.py` | Main pipeline: loads tuned params from JSON (or defaults), contains all four CV functions (`run_cv1`, `run_cv2`, `run_cv0`, `run_cross_location`), the `run_all_cv_schemes` wrapper, and CSV output saving. |
| `environment_ml.yml` | Conda environment specification (LightGBM, scikit-learn, pandas, numpy, scipy) |
| `README_env_ml.txt` | Environment setup instructions |

## Usage

```bash
cd models/LightGBM

# Step 1: Set up conda environment (once)
conda env create --file environment_ml.yml --prefix ./env-ML
conda activate ./env-ML

# Step 2: Prepare input data (once — produces all three pickles)
python Prepare_input_data.py

# Step 3: Tune hyperparameters (once per input type, produces tuned_params.json)
python tune_LightGBM.py 2>&1 | tee tune_LightGBM.out

# Step 4: Run the CV pipeline (uses tuned params if available)
python run_LightGBM.py 2>&1 | tee run_LightGBM.out
```

### Running different input types

Tuning and CV must be run **separately for each genetic feature set** (all PCs, 25 PCs, kinship). Before running steps 3 and 4, manually change the pickle path in **both** `tune_LightGBM.py` and `run_LightGBM.py`:

```python
model_input = pd.read_pickle('model_inputs/model_input.pkl')          # all 551 PCs (default)
model_input = pd.read_pickle('model_inputs/model_input_25pc.pkl')     # 25 PCs
model_input = pd.read_pickle('model_inputs/model_input_kinship.pkl')  # kinship matrix
```

**Important:** Each input type must be run from a separate working directory, because output file names (e.g. `tuned_params.json`, `cv_results_CV1_DTF_blue_LightGBM.csv`) are the same regardless of input type and would overwrite each other. Copy the scripts to separate folders (e.g. `LightGBM_allPCs/`, `LightGBM_25PCs/`, `LightGBM_kinship/`) or move outputs before running the next input type.

## Input data

### Environments

- **Locations (2):** `AUS` (Kununurra, Australia), `PAK` (Faisalabad, Pakistan)
- **Location-years (6):** `AUS_2017`, `AUS_2018`, `AUS_2019`, `PAK_2019`, `PAK_2020`, `PAK_2021`
- Not all genotypes appear in all location-years; not all traits are observed for every genotype × location-year combination (NAs present)

### Raw data files

- **Phenotype file**: `../../data/AUSPAK_phenotypes_means_BLUEs.csv`
  - Columns: `location_year`, `location`, `year`, `sample.id`, then pairs of `{trait}_mean` and `{trait}_blue` columns for each trait, plus `SdW_mean`, `SdW_z_mean`, `SdW_z_blue`
  - Only the `_blue` (BLUE) columns are used as traits; `_mean` columns are dropped by `Prepare_input_data.py`
  - Metadata columns: `location_year` (e.g. `AUS_2018`), `location` (e.g. `AUS`, `PAK`), `year`, `sample.id` (genotype identifier)
- **PCA file**: `../../data/AUSPAK_PCs_all.csv` — columns: `sample.id`, PC1-PC551
- **Kinship file**: `../../data/kinship_matrix_VanRaden_auspak_maxmissing20.csv` — VanRaden kinship matrix, sample IDs as row index and column names

### Prepared pickle files (produced by `Prepare_input_data.py`)

Each pickle contains one row per genotype × location-year observation with the following column groups:
- **Genetic features**: either `PC1`–`PC551` (all PCs), `PC1`–`PC25` (25 PCs), or `K_*` columns (kinship) — mutually exclusive per pickle
- **Trait columns**: `DTF_blue`, `DTH_blue`, `PtHt_blue`, `PcleLng_blue`, `SdLen_blue`, `TGW_blue`, `SdW_z_blue` (may contain NAs)
- **Metadata**: `sample.id`, `location_year`, `location` (retained for CV splitting/evaluation, never used as features)
- **One-hot encoded environment features**: `location_*` (e.g. `location_AUS`, `location_PAK`) and `location_year_*` (e.g. `location_year_AUS_2018`) — produced by `one_hot_encode()`, used as features by `get_feature_columns()`

Pickle files:
  - `model_inputs/model_input.pkl` — all 551 PCs
  - `model_inputs/model_input_25pc.pkl` — first 25 PCs
  - `model_inputs/model_input_kinship.pkl` — kinship columns (prefixed `K_`)

### Other input files

- **Tuned parameters**: `tuned_params.json` — per-trait LightGBM hyperparameters (produced by `tune_LightGBM.py`)

## Output files

Per-scheme, per-trait:
- `cv_results_{SCHEME}_{TRAIT}_LightGBM.csv` — accuracy metrics per location-year
- `predictions_{SCHEME}_{TRAIT}_LightGBM.csv` — individual predictions

Combined all-schemes per-trait:
- `cv_results_{TRAIT}_LightGBM_all_schemes.csv` — all CV results
- `cv_summary_{TRAIT}_LightGBM_all_schemes.csv` — summary statistics (mean/SD/min/max)
- `predictions_{TRAIT}_LightGBM_all_schemes.csv` — all predictions

## Functions shared with BayesC/RKHS/GBLUP

The following functions in `LightGBM_utils.py` are identical in logic to their counterparts in the R model pipelines:

- `VALID_TRAITS`, `LOWER_IS_BETTER_TRAITS` — trait constants
- `apply_location_year_scaling()` — z-score standardisation by location-year
- `calculate_ndcg()` — NDCG@k metric
- `evaluate_predictions()` — Pearson, Spearman, NDCG@10 evaluation
- `summarise_cv_results()` — mean/SD/min/max summary per location

## LightGBM-specific utilities

The following in `LightGBM_utils.py` are specific to the LightGBM pipeline:

- `DEFAULT_PARAMS` — pre-tuning default hyperparameters (`max_depth: 3, learning_rate: 0.05, n_estimators: 500`)
- `fetch_model_params()` — looks up per-trait hyperparameters from a tuned params dict; falls back to `DEFAULT_PARAMS` if the dict is flat or `None`
