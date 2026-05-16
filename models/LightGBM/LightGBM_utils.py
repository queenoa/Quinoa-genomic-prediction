################################################################################
# LightGBM_utils.py — Shared utility functions for LightGBM genomic prediction
#
# Sourced by:
#   run_LightGBM.py
#
# Contains:
#   - Global constants (traits, lower-is-better)
#   - Evaluation metrics (Pearson, Spearman, NDCG@10)
#   - Z-score standardisation by location-year
#   - Per-location-year evaluation
#   - CV results summarisation
#
# Functions shared with BayesC/RKHS/GBLUP (identical logic):
#   VALID_TRAITS, LOWER_IS_BETTER_TRAITS
#   apply_location_year_scaling()
#   calculate_ndcg()
#   evaluate_predictions()
#   summarise_cv_results()
################################################################################

import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import ndcg_score

# ── Global constants ─────────────────────────────────────────────────────────

VALID_TRAITS = ['DTF', 'DTH', 'PtHt', 'PcleLng',
                'SdLen', 'TGW', 'SdW_z']

LOWER_IS_BETTER_TRAITS = ['DTF', 'DTH', 'PtHt']

# ── Default model parameters ────────────────────────────────────────────────

DEFAULT_PARAMS = {'max_depth': 3, 'learning_rate': 0.05, 'n_estimators': 500}


def fetch_model_params(model_params, trait):
    """Resolve model parameters for a given trait.

    Accepts either a flat dict (used for all traits) or a dict keyed by
    trait name containing per-trait parameter dicts.
    """
    if model_params is None:
        return DEFAULT_PARAMS
    if trait in model_params and isinstance(model_params[trait], dict):
        return model_params[trait]
    return model_params


# ── NDCG@k calculation ──────────────────────────────────────────────────────

def calculate_ndcg(y_true, y_pred, k=10, lower_is_better=False):
    """Calculate NDCG@k for genomic prediction ranking quality."""
    k = min(k, len(y_true))

    if lower_is_better:
        y_true_adj = -y_true
        y_pred_adj = -y_pred
    else:
        y_true_adj = y_true
        y_pred_adj = y_pred

    # Transform to non-negative values for NDCG calculation
    min_val = min(y_true_adj.min(), y_pred_adj.min())
    if min_val < 0:
        y_true_pos = y_true_adj - min_val + 1e-6
        y_pred_pos = y_pred_adj - min_val + 1e-6
    else:
        y_true_pos = y_true_adj
        y_pred_pos = y_pred_adj

    # Reshape for sklearn (expects 2D)
    y_true_2d = y_true_pos.reshape(1, -1)
    y_pred_2d = y_pred_pos.reshape(1, -1)

    return ndcg_score(y_true_2d, y_pred_2d, k=k)

# ── Evaluate predictions ─────────────────────────────────────────────────────

def evaluate_predictions(y_true, y_pred, trait_name=None):
    """Evaluate predictions with Pearson, Spearman, and NDCG@10."""
    if len(y_true) < 2:
        return {'pearson': np.nan, 'spearman': np.nan, 'ndcg_at_10': np.nan}

    lower_is_better = trait_name in LOWER_IS_BETTER_TRAITS if trait_name else False

    try:
        rho, _ = spearmanr(y_true, y_pred)
        return {
            'pearson': np.corrcoef(y_true, y_pred)[0, 1],
            'spearman': rho,
            'ndcg_at_10': calculate_ndcg(y_true, y_pred, k=10,
                                          lower_is_better=lower_is_better)
        }
    except Exception as e:
        print(f"    Warning: Error calculating metrics: {str(e)}")
        return {'pearson': np.nan, 'spearman': np.nan, 'ndcg_at_10': np.nan}

# ── Z-score standardisation by location-year ─────────────────────────────────

def apply_location_year_scaling(data, traits):
    """Apply z-score transformation by location_year to all traits."""
    data_scaled = data.copy()

    for trait in traits:
        trait_data = data_scaled[data_scaled[trait].notna()].copy()
        if len(trait_data) < 50:
            continue

        print(f"Applying z-score transformation by location_year for {trait}...")

        for ly in trait_data['location_year'].unique():
            ly_mask = (data_scaled['location_year'] == ly) & (data_scaled[trait].notna())

            if ly_mask.sum() > 1:
                vals = data_scaled.loc[ly_mask, trait]
                ly_mean = vals.mean()
                ly_sd = vals.std()
                if ly_sd == 0 or np.isnan(ly_sd):
                    ly_sd = 1
                data_scaled.loc[ly_mask, trait] = (vals - ly_mean) / ly_sd
                print(f"  {ly}: n={ly_mask.sum()}")

    return data_scaled

# ── Evaluate per location-year ───────────────────────────────────────────────

def evaluate_per_location_year(predictions_df, trait_name, min_genotypes=10):
    """
    Evaluate predictions per location_year.

    Parameters
    ----------
    predictions_df : DataFrame
        Must have columns: sample.id, location_year, location, observed, predicted
    trait_name : str
        Trait name (for NDCG sign-flip)
    min_genotypes : int
        Minimum genotypes per location_year to compute metrics

    Returns
    -------
    DataFrame with columns: trait, location_year, location, pearson, spearman,
        ndcg_at_10, n_test_genotypes
    """
    results = []

    for ly in predictions_df['location_year'].unique():
        ly_data = predictions_df[predictions_df['location_year'] == ly]
        n_geno = ly_data['sample.id'].nunique()

        if n_geno < min_genotypes:
            print(f"      Skipping {ly} - only {n_geno} test genotypes "
                  f"(minimum: {min_genotypes})")
            continue

        y_true = ly_data['observed'].values
        y_pred = ly_data['predicted'].values

        eval_res = evaluate_predictions(y_true, y_pred, trait_name=trait_name)

        if not np.isnan(eval_res['pearson']):
            loc = ly_data['location'].iloc[0]
            results.append({
                'trait': trait_name,
                'location_year': ly,
                'location': loc,
                'pearson': eval_res['pearson'],
                'spearman': eval_res['spearman'],
                'ndcg_at_10': eval_res['ndcg_at_10'],
                'n_test_genotypes': n_geno
            })

            print(f"      {ly} - r: {eval_res['pearson']:.3f} "
                  f"| rho: {eval_res['spearman']:.3f} "
                  f"| NDCG@10: {eval_res['ndcg_at_10']:.3f} "
                  f"| N: {n_geno}")

    return pd.DataFrame(results)

# ── Summarise CV results ─────────────────────────────────────────────────────

def summarise_cv_results(cv_results, scheme_name):
    """Compute mean/SD/min/max summary per trait and location."""
    if len(cv_results) == 0:
        print(f"No successful cross-validation results for {scheme_name}")
        return pd.DataFrame()

    metric_cols = ['pearson', 'spearman', 'ndcg_at_10']

    summary = cv_results.groupby(['trait', 'location']).agg(
        **{f'{col}_{stat}': (col, stat)
           for col in metric_cols
           for stat in ['mean', 'std', 'min', 'max']},
        n_evaluations=('pearson', 'count'),
        mean_n_test_genotypes=('n_test_genotypes', 'mean')
    ).reset_index().round(4)

    print(f"\n=== {scheme_name} Cross-Validation Summary ===")
    print(summary[['trait', 'location', 'pearson_mean', 'pearson_std',
                    'spearman_mean', 'spearman_std',
                    'ndcg_at_10_mean', 'ndcg_at_10_std', 'n_evaluations']])

    return summary


print("LightGBM_utils.py loaded successfully.")
