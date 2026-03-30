# tune_LightGBM.py — Hyperparameter tuning for LightGBM genomic prediction
#
# Runs RandomizedSearchCV once per trait on the full dataset with 5-fold
# GroupKFold (grouped by genotype) to identify reasonable hyperparameters.
# Best parameters per trait are saved to tuned_params.json for use by
# run_LightGBM.py.
#
# Usage:
#   python tune_LightGBM.py

import json
import pandas as pd
import numpy as np
from lightgbm import LGBMRegressor
from sklearn.model_selection import RandomizedSearchCV, GroupKFold
from run_LightGBM import get_feature_columns
from LightGBM_utils import VALID_TRAITS, apply_location_year_scaling


PARAM_DISTRIBUTIONS = {
    'n_estimators': [100, 200, 500, 1000],
    'max_depth': [2, 3, 4, 5, 6, -1],
    'learning_rate': [0.01, 0.03, 0.05, 0.1, 0.2],
    'num_leaves': [15, 31, 63, 127],
    'min_child_samples': [5, 10, 20, 50],
    'subsample': [0.6, 0.7, 0.8, 0.9, 1.0],
    'colsample_bytree': [0.6, 0.7, 0.8, 0.9, 1.0],
    'reg_alpha': [0, 0.01, 0.1, 1.0],
    'reg_lambda': [0, 0.01, 0.1, 1.0],
}

N_ITER = 50
CV_FOLDS = 5
SCORING = 'neg_mean_squared_error'
RANDOM_STATE = 42
OUTPUT_FILE = 'tuned_params.json'


def tune_trait(data, trait, feature_columns):
    """Run RandomizedSearchCV for a single trait."""
    trait_data = data[data[trait].notna()].copy()

    if len(trait_data) < 100:
        print(f"  Skipping {trait} - only {len(trait_data)} observations")
        return None, None

    X = trait_data[feature_columns]
    y = trait_data[trait]
    groups = trait_data['sample.id']

    cv = GroupKFold(n_splits=CV_FOLDS)

    model = LGBMRegressor(random_state=RANDOM_STATE, verbosity=-1, n_jobs=1)

    search = RandomizedSearchCV(
        model,
        param_distributions=PARAM_DISTRIBUTIONS,
        n_iter=N_ITER,
        cv=cv,
        scoring=SCORING,
        random_state=RANDOM_STATE,
        n_jobs=-1,
        verbose=0,
    )

    search.fit(X, y, groups=groups)

    best = search.best_params_
    best_score = search.best_score_

    print(f"  Best neg-MSE: {best_score:.4f}")
    print(f"  Best params: {best}")

    return best, best_score


def convert_numpy_types(obj):
    """Convert numpy types to Python native types for JSON serialization."""
    if isinstance(obj, dict):
        return {k: convert_numpy_types(v) for k, v in obj.items()}
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj


if __name__ == '__main__':

    # Load and prepare data (same as run_LightGBM.py)
    model_input = pd.read_pickle('model_inputs/model_input.pkl')
    traits = VALID_TRAITS

    print("Applying z-score scaling...")
    model_input = apply_location_year_scaling(model_input, traits)

    feature_columns = get_feature_columns(model_input)
    print(f"\nFeatures: {len(feature_columns)}")
    print(f"Observations: {len(model_input)}")
    print(f"Genotypes: {model_input['sample.id'].nunique()}")
    print(f"Tuning: RandomizedSearchCV with {N_ITER} iterations, "
          f"{CV_FOLDS}-fold GroupKFold")
    print(f"Scoring: {SCORING}\n")

    tuned_params = {}
    scores = {}

    for trait in traits:
        print(f"\n--- Tuning: {trait} ---")
        best, score = tune_trait(model_input, trait, feature_columns)
        if best is not None:
            tuned_params[trait] = convert_numpy_types(best)
            scores[trait] = score

    # Save
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(tuned_params, f, indent=2)

    # Summary
    print(f"\n{'='*60}")
    print(f"Tuned parameters saved to {OUTPUT_FILE}")
    print(f"Traits tuned: {len(tuned_params)} / {len(traits)}")
    print(f"{'='*60}")
    print(f"\n{'Trait':<16} {'neg-MSE':>10}")
    print(f"{'-'*16} {'-'*10}")
    for trait in traits:
        if trait in scores:
            print(f"{trait:<16} {scores[trait]:>10.4f}")
