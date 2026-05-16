# LightGBM cross-validation with z-score transformation
# Extended with multiple CV schemes:
#   CV1: New genotypes in known location-years
#   CV2: Sparse testing (known genotypes, incomplete location-years)
#   CV0: Leave-one-location-year-out
#   Cross-location transferability: train on one location -> predict others
#
# Structure mirrors GBLUP.R — all CV schemes in a single script.

import json
import os
import pandas as pd
import numpy as np
from lightgbm import LGBMRegressor
from sklearn.model_selection import GroupKFold, StratifiedKFold
from LightGBM_utils import (
    VALID_TRAITS, LOWER_IS_BETTER_TRAITS, DEFAULT_PARAMS,
    apply_location_year_scaling, evaluate_predictions,
    evaluate_per_location_year, summarise_cv_results,
    fetch_model_params
)


# ============================================================================
# Feature column helpers
# ============================================================================

def get_feature_columns(data, include_location=True, include_location_year=True):
    """Get feature column names based on what should be included.

    Genetic features are kinship matrix columns (prefix 'K_').
    """
    cols = [c for c in data.columns if c.startswith('K_')]
    if include_location:
        cols += [c for c in data.columns if c.startswith('location_') and
                 not c.startswith('location_year')]
    if include_location_year:
        cols += [c for c in data.columns if c.startswith('location_year_')]
    return cols


# ============================================================================
# CV1: New genotypes in known location-years
# Test genotypes have their trait values masked across all location-years.
# ============================================================================

def run_cv1(data, traits, model_class=LGBMRegressor, model_params=None,
            k_folds=5, n_iterations=15, min_genotypes=10):

    print("\n==========================================")
    print("CV1: New genotypes in known location-years")
    print("==========================================\n")

    if model_params is None:
        model_params = DEFAULT_PARAMS

    feature_columns = get_feature_columns(data)
    genotypes = data['sample.id'].unique()
    n_genotypes = len(genotypes)

    print(f"Data: {n_genotypes} genotypes, "
          f"{data['location_year'].nunique()} location-years, "
          f"{len(data)} observations")
    print(f"Settings: {k_folds} folds, {n_iterations} iterations, "
          f"min genotypes per location-year: {min_genotypes}")
    print(f"Features: {len(feature_columns)} (kinship + encoded variables)\n")

    cv_results = []
    cv_predictions = []

    for iter_num in range(1, n_iterations + 1):
        current_seed = 1000 + iter_num
        print(f"Iteration {iter_num} of {n_iterations} (seed: {current_seed})")

        gkf = GroupKFold(n_splits=k_folds, shuffle=True, random_state=current_seed)
        fold_splits = list(gkf.split(data, groups=data['sample.id']))

        for fold_idx, (train_idx, test_idx) in enumerate(fold_splits, 1):
            print(f"  Fold {fold_idx} of {k_folds}")

            test_genotypes = data.iloc[test_idx]['sample.id'].unique()

            for trait in traits:
                try:
                    trait_data = data[data[trait].notna()].copy()
                    if len(trait_data) < 50:
                        continue

                    # Split by genotype
                    train_data = trait_data[~trait_data['sample.id'].isin(test_genotypes)]
                    test_data = trait_data[trait_data['sample.id'].isin(test_genotypes)]

                    if len(train_data) < 50:
                        continue

                    # Train model
                    X_train = train_data[feature_columns]
                    y_train = train_data[trait]
                    params = fetch_model_params(model_params, trait)
                    model = model_class(random_state=current_seed,
                                        verbosity=-1, **params)
                    model.fit(X_train, y_train)

                    # Predict and evaluate per location_year
                    X_test = test_data[feature_columns]
                    test_data = test_data.copy()
                    test_data['predicted'] = model.predict(X_test)

                    # Average by genotype within each location_year
                    geno_avg = test_data.groupby(
                        ['sample.id', 'location_year', 'location']
                    ).agg({trait: 'mean', 'predicted': 'mean'}).reset_index()
                    geno_avg = geno_avg.rename(columns={trait: 'observed'})

                    # Collect predictions
                    preds = geno_avg[['sample.id', 'location_year', 'location', 'observed', 'predicted']].copy()
                    preds['trait'] = trait
                    preds['iteration'] = iter_num
                    preds['fold'] = fold_idx
                    preds['seed'] = current_seed
                    preds['cv_scheme'] = 'CV1'
                    cv_predictions.append(preds)

                    # Evaluate per location_year
                    print(f"    Trait: {trait}")
                    ly_metrics = evaluate_per_location_year(
                        geno_avg, trait, min_genotypes=min_genotypes
                    )

                    if len(ly_metrics) > 0:
                        ly_metrics['iteration'] = iter_num
                        ly_metrics['fold'] = fold_idx
                        ly_metrics['seed'] = current_seed
                        ly_metrics['cv_scheme'] = 'CV1'
                        cv_results.append(ly_metrics)

                except Exception as e:
                    print(f"    Error in trait {trait}: {str(e)}")
                    continue

    cv_results = pd.concat(cv_results, ignore_index=True) if cv_results else pd.DataFrame()
    cv_predictions = pd.concat(cv_predictions, ignore_index=True) if cv_predictions else pd.DataFrame()
    summary = summarise_cv_results(cv_results, "CV1")

    return {
        'results': cv_results,
        'predictions': cv_predictions,
        'summary': summary,
        'cv_scheme': 'CV1'
    }


# ============================================================================
# CV2: Sparse testing (5-fold stratified by location-year)
# Observed cells are assigned to 5 folds within each location-year.
# Masking is at the observation level, not the genotype level.
# ============================================================================

def run_cv2(data, traits, model_class=LGBMRegressor, model_params=None,
            k_folds=5, n_iterations=15, min_genotypes=10):

    print("\n==================================================")
    print("CV2: Sparse testing (5-fold stratified)")
    print("==================================================\n")

    if model_params is None:
        model_params = DEFAULT_PARAMS

    feature_columns = get_feature_columns(data)

    print(f"Data: {data['sample.id'].nunique()} genotypes, "
          f"{data['location_year'].nunique()} location-years, "
          f"{len(data)} observations")
    print(f"Settings: {k_folds} folds, {n_iterations} iterations, "
          f"min genotypes per location-year: {min_genotypes}")
    print(f"Features: {len(feature_columns)} (kinship + encoded variables)\n")

    cv_results = []
    cv_predictions = []

    for iter_num in range(1, n_iterations + 1):
        current_seed = 2000 + iter_num

        print(f"Iteration {iter_num} of {n_iterations} (seed: {current_seed})")

        for trait in traits:
            try:
                trait_data = data[data[trait].notna()].copy()
                n_observed = len(trait_data)

                if n_observed < 100:
                    print(f"  Too few observations for {trait}: {n_observed}")
                    continue

                skf = StratifiedKFold(n_splits=k_folds, shuffle=True,
                                      random_state=current_seed)

                print(f"  Trait: {trait}")

                for fold, (train_idx, test_idx) in enumerate(
                    skf.split(trait_data, trait_data['location_year']), 1
                ):
                    print(f"    Fold {fold} of {k_folds}")

                    train_set = trait_data.iloc[train_idx]
                    test_set = trait_data.iloc[test_idx].copy()

                    if len(train_set) < 50:
                        print(f"      Too few training observations: {len(train_set)}")
                        continue

                    # Train model
                    X_train = train_set[feature_columns]
                    y_train = train_set[trait]
                    params = fetch_model_params(model_params, trait)
                    model = model_class(random_state=current_seed,
                                        verbosity=-1, **params)
                    model.fit(X_train, y_train)

                    # Predict test cells
                    X_test = test_set[feature_columns]
                    test_set['predicted'] = model.predict(X_test)
                    test_set = test_set.rename(columns={trait: 'observed'})

                    # Collect predictions
                    preds = test_set[['sample.id', 'location_year', 'location', 'observed', 'predicted']].copy()
                    preds['trait'] = trait
                    preds['iteration'] = iter_num
                    preds['fold'] = fold
                    preds['seed'] = current_seed
                    preds['cv_scheme'] = 'CV2'
                    cv_predictions.append(preds)

                    # Evaluate per location_year
                    pred_df = test_set[['sample.id', 'location_year', 'location',
                                        'observed', 'predicted']].copy()
                    ly_metrics = evaluate_per_location_year(
                        pred_df, trait, min_genotypes=min_genotypes
                    )

                    if len(ly_metrics) > 0:
                        ly_metrics['iteration'] = iter_num
                        ly_metrics['fold'] = fold
                        ly_metrics['seed'] = current_seed
                        ly_metrics['cv_scheme'] = 'CV2'
                        cv_results.append(ly_metrics)

            except Exception as e:
                print(f"  Error in trait {trait}: {str(e)}")
                continue

    cv_results = pd.concat(cv_results, ignore_index=True) if cv_results else pd.DataFrame()
    cv_predictions = pd.concat(cv_predictions, ignore_index=True) if cv_predictions else pd.DataFrame()
    summary = summarise_cv_results(cv_results, "CV2")

    return {
        'results': cv_results,
        'predictions': cv_predictions,
        'summary': summary,
        'cv_scheme': 'CV2'
    }


# ============================================================================
# CV0: Leave-one-location-year-out
# All rows kept; trait masked for the held-out location-year.
# Deterministic — no random iterations.
# ============================================================================

def run_cv0(data, traits, model_class=LGBMRegressor, model_params=None,
            min_genotypes=10):

    print("\n=============================================")
    print("CV0: Leave-one-location-year-out")
    print("=============================================\n")

    if model_params is None:
        model_params = DEFAULT_PARAMS

    feature_columns = get_feature_columns(data)
    location_years = data['location_year'].unique()

    print(f"Data: {data['sample.id'].nunique()} genotypes, "
          f"{len(location_years)} location-years, "
          f"{len(data)} observations")
    print(f"Min genotypes: {min_genotypes}")
    print(f"Features: {len(feature_columns)} (kinship + encoded variables)\n")

    cv_results = []
    cv_predictions = []

    for held_out_ly in location_years:
        print(f"Holding out: {held_out_ly}")

        for trait in traits:
            try:
                trait_data = data[data[trait].notna()].copy()

                # Test set: observations in held-out location_year
                test_set = trait_data[trait_data['location_year'] == held_out_ly].copy()
                # Training set: all other location_years
                train_set = trait_data[trait_data['location_year'] != held_out_ly]

                test_genos = test_set['sample.id'].unique()
                if len(test_genos) < min_genotypes:
                    print(f"  Skipping {trait} - only {len(test_genos)} "
                          f"genotypes in held-out location-year")
                    continue

                if len(train_set) < 50:
                    continue

                # Train model
                X_train = train_set[feature_columns]
                y_train = train_set[trait]
                params = fetch_model_params(model_params, trait)
                model = model_class(random_state=42, verbosity=-1, **params)
                model.fit(X_train, y_train)

                # Predict held-out location_year
                X_test = test_set[feature_columns]
                test_set['predicted'] = model.predict(X_test)

                # Average by genotype within the held-out location_year
                geno_avg = test_set.groupby(
                    ['sample.id', 'location_year', 'location']
                ).agg({trait: 'mean', 'predicted': 'mean'}).reset_index()
                geno_avg = geno_avg.rename(columns={trait: 'observed'})

                # Collect predictions
                preds = geno_avg[['sample.id', 'location_year', 'location', 'observed', 'predicted']].copy()
                preds['trait'] = trait
                preds['cv_scheme'] = 'CV0'
                cv_predictions.append(preds)

                # Evaluate
                print(f"  Trait: {trait}")
                ly_metrics = evaluate_per_location_year(
                    geno_avg, trait, min_genotypes=min_genotypes
                )

                if len(ly_metrics) > 0:
                    ly_metrics['iteration'] = np.nan
                    ly_metrics['fold'] = np.nan
                    ly_metrics['seed'] = np.nan
                    ly_metrics['cv_scheme'] = 'CV0'
                    cv_results.append(ly_metrics)

            except Exception as e:
                print(f"  Error in trait {trait} for {held_out_ly}: {str(e)}")
                continue

    cv_results = pd.concat(cv_results, ignore_index=True) if cv_results else pd.DataFrame()
    cv_predictions = pd.concat(cv_predictions, ignore_index=True) if cv_predictions else pd.DataFrame()
    summary = summarise_cv_results(cv_results, "CV0")

    return {
        'results': cv_results,
        'predictions': cv_predictions,
        'summary': summary,
        'cv_scheme': 'CV0'
    }


# ============================================================================
# Cross-location transferability
# Train on one location only, predict all other locations.
# Uses kinship features only (no location/location_year encoding) —
# prediction comes purely from genetic relationships, matching GBLUP's
# trait ~ 1 + G approach.
# ============================================================================

def run_cross_location(data, traits, model_class=LGBMRegressor, model_params=None,
                        min_genotypes=10):

    print("\n=============================================")
    print("Cross-location transferability")
    print("(train on one location -> predict others)")
    print("=============================================\n")

    if model_params is None:
        model_params = DEFAULT_PARAMS

    # Kinship only — no location/location_year encoding
    genetic_columns = get_feature_columns(data, include_location=False,
                                           include_location_year=False)
    locations = data['location'].unique()

    print(f"Data: {data['sample.id'].nunique()} genotypes, "
          f"{data['location_year'].nunique()} location-years, "
          f"{len(data)} observations")
    print(f"Locations: {', '.join(locations)}")
    print(f"Features: {len(genetic_columns)} (kinship only)")
    print(f"Min genotypes: {min_genotypes}\n")

    cv_results = []
    cv_predictions = []

    for train_loc in locations:
        predict_locs = [loc for loc in locations if loc != train_loc]
        print(f"Training on: {train_loc} -> Predicting: {', '.join(predict_locs)}")

        for trait in traits:
            try:
                trait_data = data[data[trait].notna()].copy()
                train_data = trait_data[trait_data['location'] == train_loc]

                n_train_obs = len(train_data)
                if n_train_obs < 50:
                    print(f"  Skipping {trait} - only {n_train_obs} training observations")
                    continue

                # Train model on source location only
                X_train = train_data[genetic_columns]
                y_train = train_data[trait]
                params = fetch_model_params(model_params, trait)
                model = model_class(random_state=42, verbosity=-1, **params)
                model.fit(X_train, y_train)

                cv_label = f"CrossLoc_{train_loc}->{'+'.join(predict_locs)}"

                # Predict each target location_year
                target_data = trait_data[trait_data['location'].isin(predict_locs)]
                target_lys = target_data['location_year'].unique()

                print(f"  Trait: {trait}")

                for ly in target_lys:
                    ly_data = target_data[target_data['location_year'] == ly].copy()

                    if len(ly_data) < 3:
                        continue

                    # Predict
                    X_test = ly_data[genetic_columns]
                    ly_data['predicted'] = model.predict(X_test)

                    # Average by genotype
                    geno_avg = ly_data.groupby(
                        ['sample.id', 'location_year', 'location']
                    ).agg({trait: 'mean', 'predicted': 'mean'}).reset_index()
                    geno_avg = geno_avg.rename(columns={trait: 'observed'})

                    n_geno = len(geno_avg)

                    # Collect predictions
                    preds = geno_avg[['sample.id', 'location_year', 'location', 'observed', 'predicted']].copy()
                    preds['trait'] = trait
                    preds['train_location'] = train_loc
                    preds['cv_scheme'] = cv_label
                    cv_predictions.append(preds)

                    if n_geno < min_genotypes:
                        print(f"      Skipping {ly} - only {n_geno} "
                              f"overlapping genotypes (min: {min_genotypes})")
                        continue

                    eval_res = evaluate_predictions(
                        geno_avg['observed'].values,
                        geno_avg['predicted'].values,
                        trait_name=trait
                    )

                    if not np.isnan(eval_res['pearson']):
                        cv_results.append({
                            'iteration': np.nan,
                            'fold': np.nan,
                            'trait': trait,
                            'location_year': ly,
                            'location': geno_avg['location'].iloc[0],
                            'pearson': eval_res['pearson'],
                            'spearman': eval_res['spearman'],
                            'ndcg_at_10': eval_res['ndcg_at_10'],
                            'seed': np.nan,
                            'n_test_genotypes': n_geno,
                            'cv_scheme': cv_label
                        })

                        print(f"      {ly} - r: {eval_res['pearson']:.3f} "
                              f"| rho: {eval_res['spearman']:.3f} "
                              f"| NDCG@10: {eval_res['ndcg_at_10']:.3f} "
                              f"| N: {n_geno}")

            except Exception as e:
                print(f"  Error in trait {trait}: {str(e)}")
                continue

    cv_results = pd.DataFrame(cv_results)
    cv_predictions = pd.concat(cv_predictions, ignore_index=True) if cv_predictions else pd.DataFrame()

    # Summary grouped by direction
    if len(cv_results) > 0:
        summary = cv_results.groupby(['trait', 'cv_scheme']).agg(
            pearson_mean=('pearson', 'mean'),
            pearson_std=('pearson', 'std'),
            pearson_min=('pearson', 'min'),
            pearson_max=('pearson', 'max'),
            spearman_mean=('spearman', 'mean'),
            spearman_std=('spearman', 'std'),
            spearman_min=('spearman', 'min'),
            spearman_max=('spearman', 'max'),
            ndcg_at_10_mean=('ndcg_at_10', 'mean'),
            ndcg_at_10_std=('ndcg_at_10', 'std'),
            ndcg_at_10_min=('ndcg_at_10', 'min'),
            ndcg_at_10_max=('ndcg_at_10', 'max'),
            n_target_location_years=('pearson', 'count'),
            mean_n_test_genotypes=('n_test_genotypes', 'mean')
        ).reset_index().round(4)

        print("\n=== Cross-Location Transferability Summary ===")
        print(summary[['trait', 'cv_scheme', 'pearson_mean', 'pearson_std',
                        'spearman_mean', 'spearman_std',
                        'ndcg_at_10_mean', 'ndcg_at_10_std',
                        'n_target_location_years']])
    else:
        print("No successful cross-location results")
        summary = pd.DataFrame()

    return {
        'results': cv_results,
        'predictions': cv_predictions,
        'summary': summary,
        'cv_scheme': 'CrossLocation'
    }


# ============================================================================
# WRAPPER: Run all CV schemes and combine results
# ============================================================================

def run_all_cv_schemes(data, traits, model_class=LGBMRegressor, model_params=None,
                        k_folds=5, n_iterations=15, min_genotypes=10):

    print("==========================================================")
    print("Running all cross-validation schemes for LightGBM")
    print("==========================================================\n")

    results = {}

    results['cv1'] = run_cv1(data, traits, model_class, model_params,
                              k_folds=k_folds, n_iterations=n_iterations,
                              min_genotypes=min_genotypes)

    results['cv2'] = run_cv2(data, traits, model_class, model_params,
                              k_folds=k_folds, n_iterations=n_iterations,
                              min_genotypes=min_genotypes)

    results['cv0'] = run_cv0(data, traits, model_class, model_params,
                              min_genotypes=min_genotypes)

    results['cross_loc'] = run_cross_location(data, traits, model_class,
                                               model_params,
                                               min_genotypes=min_genotypes)

    # Combine all results
    all_results_list = [results[k]['results'] for k in results
                        if len(results[k]['results']) > 0]
    all_results = pd.concat(all_results_list, ignore_index=True) if all_results_list else pd.DataFrame()

    print("\n==========================================================")
    print("All CV schemes complete")
    print(f"Total evaluations: {len(all_results)}")
    if len(all_results) > 0:
        print("Breakdown:")
        print(all_results['cv_scheme'].value_counts())
    print("==========================================================")

    results['all_results'] = all_results
    return results


# ============================================================================
# RUN PIPELINE
# ============================================================================

if __name__ == '__main__':

    # Load tuned hyperparameters if available, otherwise use defaults
    params_file = 'tuned_params.json'
    if os.path.exists(params_file):
        with open(params_file) as f:
            model_params = json.load(f)
        print(f"Loaded tuned hyperparameters from {params_file}")
        for trait_name, trait_params in model_params.items():
            print(f"  {trait_name}: {trait_params}")
    else:
        model_params = DEFAULT_PARAMS
        print(f"Using default hyperparameters (no {params_file} found)")
        print(f"  {model_params}")

    # Load model input data
    model_input = pd.read_pickle('model_inputs/model_input.pkl')

    # Apply z-score scaling upfront
    traits = VALID_TRAITS
    print(f"Traits: {', '.join(traits)}")
    print(f"Z-score transformation: ENABLED\n")
    model_input = apply_location_year_scaling(model_input, traits)

    # Run all CV schemes
    all_cv = run_all_cv_schemes(
        data=model_input,
        traits=traits,
        model_class=LGBMRegressor,
        model_params=model_params,
        k_folds=5,
        n_iterations=15,
        min_genotypes=10
    )

    # ========================================================================
    # SAVE OUTPUTS
    # ========================================================================

    scheme_labels = {
        'cv1': 'CV1',
        'cv2': 'CV2',
        'cv0': 'CV0',
        'cross_loc': 'CrossLoc'
    }

    # --- Per-scheme, per-trait CSVs ---
    for scheme, label in scheme_labels.items():
        res = all_cv[scheme]['results']
        preds = all_cv[scheme]['predictions']

        for trait in traits:
            if len(res) > 0:
                trait_res = res[res['trait'] == trait]
                if len(trait_res) > 0:
                    outfile = f"cv_results_{label}_{trait}_LightGBM.csv"
                    trait_res.to_csv(outfile, index=False)
                    print(f"Saved: {outfile} ({len(trait_res)} rows)")

            if len(preds) > 0:
                trait_preds = preds[preds['trait'] == trait]
                if len(trait_preds) > 0:
                    outfile = f"predictions_{label}_{trait}_LightGBM.csv"
                    trait_preds.to_csv(outfile, index=False)
                    print(f"Saved: {outfile} ({len(trait_preds)} rows)")

    # --- Combined all-schemes per-trait CSVs ---
    all_preds_list = [all_cv[k]['predictions'] for k in scheme_labels
                      if len(all_cv[k]['predictions']) > 0]
    all_preds = pd.concat(all_preds_list, ignore_index=True) if all_preds_list else pd.DataFrame()

    all_summaries_list = []
    for scheme, label in scheme_labels.items():
        s = all_cv[scheme]['summary']
        if len(s) > 0:
            s = s.copy()
            s['cv_scheme'] = label
            all_summaries_list.append(s)
    all_summaries = pd.concat(all_summaries_list, ignore_index=True) if all_summaries_list else pd.DataFrame()

    for trait in traits:
        # All results across schemes
        if len(all_cv['all_results']) > 0:
            trait_all = all_cv['all_results'][all_cv['all_results']['trait'] == trait]
            if len(trait_all) > 0:
                outfile = f"cv_results_{trait}_LightGBM_all_schemes.csv"
                trait_all.to_csv(outfile, index=False)
                print(f"Saved: {outfile} ({len(trait_all)} rows)")

        # Summary statistics across schemes
        if len(all_summaries) > 0:
            trait_summary = all_summaries[all_summaries['trait'] == trait]
            if len(trait_summary) > 0:
                outfile = f"cv_summary_{trait}_LightGBM_all_schemes.csv"
                trait_summary.to_csv(outfile, index=False)
                print(f"Saved: {outfile} ({len(trait_summary)} rows)")

        # All predictions across schemes
        if len(all_preds) > 0:
            trait_preds = all_preds[all_preds['trait'] == trait]
            if len(trait_preds) > 0:
                outfile = f"predictions_{trait}_LightGBM_all_schemes.csv"
                trait_preds.to_csv(outfile, index=False)
                print(f"Saved: {outfile} ({len(trait_preds)} rows)")

    print("\n=== All outputs saved ===")
