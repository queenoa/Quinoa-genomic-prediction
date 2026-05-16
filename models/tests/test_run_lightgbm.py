"""Tests for run_LightGBM.py — feature columns and CV scheme functions.

Uses a 100-genotype subset of real AUSPAK data via fixtures.
A lightweight LinearRegression stands in for LGBMRegressor to keep tests fast.

Production pipeline is kinship-only (commit c97e97a dropped PC features), so
all tests use the kinship-feature fixture.
"""

import numpy as np
import pandas as pd
import pytest
from sklearn.linear_model import LinearRegression

from run_LightGBM import (
    get_feature_columns, run_cv1, run_cv2, run_cv0, run_cross_location,
    run_all_cv_schemes,
)
from LightGBM_utils import VALID_TRAITS, DEFAULT_PARAMS, apply_location_year_scaling


class _QuickRegressor(LinearRegression):
    """LinearRegression wrapper that accepts (and ignores) LightGBM kwargs.

    Stores all kwargs as attributes so sklearn's get_params()/clone() work.
    """

    def __init__(self, random_state=None, verbosity=None, **kwargs):
        self.random_state = random_state
        self.verbosity = verbosity
        for k, v in kwargs.items():
            setattr(self, k, v)
        super().__init__()


# ── get_feature_columns ─────────────────────────────────────────────────────

class TestGetFeatureColumns:

    def test_kinship_features_detected(self, model_input_kinship):
        cols = get_feature_columns(model_input_kinship)
        k_cols = [c for c in cols if c.startswith('K_')]
        assert len(k_cols) == 100  # one per genotype in subset

    def test_location_included_by_default(self, model_input_kinship):
        """include_location adds location_AUS, location_PAK columns."""
        cols = get_feature_columns(model_input_kinship)
        pure_loc = [c for c in cols if c.startswith('location_')
                    and not c.startswith('location_year')]
        assert len(pure_loc) == 2  # AUS, PAK

    def test_location_year_included_by_default(self, model_input_kinship):
        cols = get_feature_columns(model_input_kinship)
        ly_cols = [c for c in cols if c.startswith('location_year_')]
        assert len(ly_cols) > 0

    def test_exclude_location_removes_pure_location(self, model_input_kinship):
        """With include_location=False, location_AUS/PAK are absent."""
        cols = get_feature_columns(model_input_kinship, include_location=False)
        pure_loc = [c for c in cols if c.startswith('location_')
                    and not c.startswith('location_year')]
        assert len(pure_loc) == 0

    def test_exclude_both_removes_all_location_features(self, model_input_kinship):
        cols = get_feature_columns(model_input_kinship, include_location=False,
                                   include_location_year=False)
        loc_cols = [c for c in cols if c.startswith('location')]
        assert len(loc_cols) == 0

    def test_genetic_only_returns_kinship(self, model_input_kinship):
        cols = get_feature_columns(model_input_kinship, include_location=False,
                                   include_location_year=False)
        assert all(c.startswith('K_') for c in cols)

    def test_pc_columns_not_detected(self, model_input_kinship):
        """get_feature_columns should not pick up PC* columns (PCs dropped from
        production in commit c97e97a)."""
        df = model_input_kinship.copy()
        df['PC1'] = 0.0
        df['PC2'] = 0.0
        cols = get_feature_columns(df)
        assert 'PC1' not in cols
        assert 'PC2' not in cols

    def test_no_trait_or_metadata_columns(self, model_input_kinship):
        cols = get_feature_columns(model_input_kinship)
        assert 'sample.id' not in cols
        assert 'DTF' not in cols
        assert 'location' not in cols
        assert 'location_year' not in cols

    def test_no_duplicate_columns(self, model_input_kinship):
        """Feature list should have no duplicates (LightGBM rejects them)."""
        cols = get_feature_columns(model_input_kinship)
        assert len(cols) == len(set(cols))


# ── Scaled fixture for CV tests ─────────────────────────────────────────────

@pytest.fixture
def scaled_kinship(model_input_kinship):
    traits = [t for t in VALID_TRAITS if t in model_input_kinship.columns]
    return apply_location_year_scaling(model_input_kinship, traits)


# ── CV1: New genotypes ──────────────────────────────────────────────────────

class TestRunCV1:

    def test_returns_expected_keys(self, scaled_kinship):
        out = run_cv1(scaled_kinship, ['TGW'],
                      model_class=_QuickRegressor,
                      model_params=DEFAULT_PARAMS,
                      k_folds=3, n_iterations=2, min_genotypes=3)
        assert set(out.keys()) == {'results', 'predictions', 'summary',
                                   'cv_scheme'}
        assert out['cv_scheme'] == 'CV1'

    def test_predictions_have_expected_columns(self, scaled_kinship):
        out = run_cv1(scaled_kinship, ['TGW'],
                      model_class=_QuickRegressor,
                      k_folds=3, n_iterations=1, min_genotypes=3)
        expected = {'sample.id', 'location_year', 'location', 'observed',
                    'predicted', 'trait', 'iteration', 'fold', 'seed',
                    'cv_scheme'}
        assert expected == set(out['predictions'].columns)

    def test_all_genotypes_with_data_predicted(self, scaled_kinship):
        """Every genotype with non-null trait values should appear in predictions."""
        trait = 'TGW'
        out = run_cv1(scaled_kinship, [trait],
                      model_class=_QuickRegressor,
                      k_folds=5, n_iterations=1, min_genotypes=3)
        predicted_genos = set(out['predictions']['sample.id'].unique())
        genos_with_data = set(
            scaled_kinship.loc[scaled_kinship[trait].notna(), 'sample.id'].unique()
        )
        assert predicted_genos == genos_with_data

    def test_seed_scheme_starts_at_1001(self, scaled_kinship):
        out = run_cv1(scaled_kinship, ['TGW'],
                      model_class=_QuickRegressor,
                      k_folds=3, n_iterations=2, min_genotypes=3)
        seeds = out['predictions']['seed'].unique()
        assert 1001 in seeds
        assert 1002 in seeds

    def test_folds_differ_across_iterations(self, scaled_kinship):
        """Different iterations should assign genotypes to different folds."""
        trait = 'TGW'
        out = run_cv1(scaled_kinship, [trait],
                      model_class=_QuickRegressor,
                      k_folds=3, n_iterations=2, min_genotypes=3)
        preds = out['predictions']
        iter1_fold1 = set(preds[(preds['iteration'] == 1) & (preds['fold'] == 1)]['sample.id'])
        iter2_fold1 = set(preds[(preds['iteration'] == 2) & (preds['fold'] == 1)]['sample.id'])
        assert iter1_fold1 != iter2_fold1, (
            "Fold 1 has identical genotypes in iterations 1 and 2 — "
            "GroupKFold may not be receiving shuffled data"
        )


# ── CV2: Sparse testing ─────────────────────────────────────────────────────

class TestRunCV2:

    def test_returns_expected_keys(self, scaled_kinship):
        out = run_cv2(scaled_kinship, ['TGW'],
                      model_class=_QuickRegressor,
                      k_folds=3, n_iterations=1, min_genotypes=3)
        assert out['cv_scheme'] == 'CV2'
        assert 'results' in out

    def test_seed_scheme_starts_at_2001(self, scaled_kinship):
        out = run_cv2(scaled_kinship, ['TGW'],
                      model_class=_QuickRegressor,
                      k_folds=3, n_iterations=2, min_genotypes=3)
        seeds = out['predictions']['seed'].unique()
        assert 2001 in seeds
        assert 2002 in seeds

    def test_observation_level_masking(self, scaled_kinship):
        """CV2 masks at observation level — total predictions ≈ non-null obs."""
        out = run_cv2(scaled_kinship, ['TGW'],
                      model_class=_QuickRegressor,
                      k_folds=3, n_iterations=1, min_genotypes=3)
        n_obs = scaled_kinship['TGW'].notna().sum()
        n_pred = len(out['predictions'])
        assert n_pred == n_obs


# ── CV0: Leave-one-location-year-out ────────────────────────────────────────

class TestRunCV0:

    def test_returns_expected_keys(self, scaled_kinship):
        out = run_cv0(scaled_kinship, ['TGW'],
                      model_class=_QuickRegressor,
                      min_genotypes=3)
        assert out['cv_scheme'] == 'CV0'
        assert 'results' in out

    def test_deterministic_no_random_iterations(self, scaled_kinship):
        """CV0 has no iteration/fold — these should be NaN."""
        out = run_cv0(scaled_kinship, ['TGW'],
                      model_class=_QuickRegressor,
                      min_genotypes=3)
        if len(out['results']) > 0:
            assert out['results']['iteration'].isna().all()
            assert out['results']['fold'].isna().all()

    def test_each_location_year_with_data_held_out(self, scaled_kinship):
        """Each location-year that has non-null trait data should appear in predictions."""
        trait = 'TGW'
        out = run_cv0(scaled_kinship, [trait],
                      model_class=_QuickRegressor,
                      min_genotypes=3)
        ly_in_preds = set(out['predictions']['location_year'].unique())
        ly_with_data = set(
            scaled_kinship.loc[scaled_kinship[trait].notna(), 'location_year'].unique()
        )
        assert ly_in_preds == ly_with_data


# ── Cross-location ──────────────────────────────────────────────────────────

class TestRunCrossLocation:

    def test_returns_expected_keys(self, scaled_kinship):
        out = run_cross_location(scaled_kinship, ['TGW'],
                                 model_class=_QuickRegressor,
                                 min_genotypes=3)
        assert out['cv_scheme'] == 'CrossLocation'

    def test_uses_kinship_features_only(self, scaled_kinship):
        """CrossLoc should use only K_ features, no location encoding."""
        cols = get_feature_columns(scaled_kinship, include_location=False,
                                   include_location_year=False)
        assert all(c.startswith('K_') for c in cols)
        assert not any(c.startswith('location_') for c in cols)
        assert not any(c.startswith('location_year_') for c in cols)

    def test_train_location_not_in_predictions(self, scaled_kinship):
        out = run_cross_location(scaled_kinship, ['TGW'],
                                 model_class=_QuickRegressor,
                                 min_genotypes=3)
        preds = out['predictions']
        if len(preds) > 0 and 'train_location' in preds.columns:
            for _, row in preds.iterrows():
                assert row['location'] != row['train_location']


# ── Evaluation counts per CV scheme ────────────────────────────────────────
# These tests verify that each scheme produces the expected number of
# per-location-year evaluations.  The CV2 fold-assignment bug (all
# observations in a location-year collapsed to one fold) would cause CV2
# to produce ~1/k_folds the evaluations of CV1.

def _n_qualifying_lys(data, trait, min_genotypes):
    """Count location-years with enough non-null observations for evaluation."""
    trait_data = data[data[trait].notna()]
    return (trait_data.groupby('location_year')['sample.id'].nunique()
            >= min_genotypes).sum()


class TestEvaluationCounts:

    def test_cv1_evaluation_count(self, scaled_kinship):
        """CV1: n_iterations × k_folds × n_qualifying_location_years per trait."""
        trait = 'TGW'
        k_folds, n_iter, min_g = 3, 2, 3
        out = run_cv1(scaled_kinship, [trait], model_class=_QuickRegressor,
                      k_folds=k_folds, n_iterations=n_iter, min_genotypes=min_g)
        n_ly = _n_qualifying_lys(scaled_kinship, trait, min_g)
        expected = n_iter * k_folds * n_ly
        assert len(out['results']) == expected

    def test_cv2_evaluation_count_matches_cv1(self, scaled_kinship):
        """CV2 should produce the same number of evaluations as CV1.

        Both use k_folds × n_iterations with per-location-year evaluation.
        The bug where np.random.choice returned a scalar (assigning all
        observations in a location-year to one fold) would make CV2 produce
        roughly 1/k_folds the evaluations of CV1.
        """
        trait = 'TGW'
        k_folds, n_iter, min_g = 3, 2, 3
        cv1 = run_cv1(scaled_kinship, [trait], model_class=_QuickRegressor,
                      k_folds=k_folds, n_iterations=n_iter, min_genotypes=min_g)
        cv2 = run_cv2(scaled_kinship, [trait], model_class=_QuickRegressor,
                      k_folds=k_folds, n_iterations=n_iter, min_genotypes=min_g)
        assert len(cv2['results']) == len(cv1['results'])

    def test_cv2_all_location_years_in_every_fold(self, scaled_kinship):
        """Each fold in CV2 should contain evaluations for all qualifying
        location-years — not just a random subset."""
        trait = 'TGW'
        k_folds, min_g = 3, 3
        out = run_cv2(scaled_kinship, [trait], model_class=_QuickRegressor,
                      k_folds=k_folds, n_iterations=1, min_genotypes=min_g)
        results = out['results']
        n_ly = _n_qualifying_lys(scaled_kinship, trait, min_g)

        for fold in range(1, k_folds + 1):
            fold_lys = results[results['fold'] == fold]['location_year'].nunique()
            assert fold_lys == n_ly, (
                f"Fold {fold} evaluated {fold_lys} location-years, expected "
                f"{n_ly}. Fold assignment may not distribute observations "
                f"across folds within each location-year."
            )

    def test_cv0_evaluation_count(self, scaled_kinship):
        """CV0: exactly one evaluation per qualifying location-year per trait."""
        trait = 'TGW'
        min_g = 3
        out = run_cv0(scaled_kinship, [trait], model_class=_QuickRegressor,
                      min_genotypes=min_g)
        n_ly = _n_qualifying_lys(scaled_kinship, trait, min_g)
        assert len(out['results']) == n_ly

    def test_crossloc_evaluation_count(self, scaled_kinship):
        """CrossLoc: for each training location, one evaluation per qualifying
        target location-year."""
        trait = 'TGW'
        min_g = 3
        out = run_cross_location(scaled_kinship, [trait],
                                 model_class=_QuickRegressor,
                                 min_genotypes=min_g)
        results = out['results']
        trait_data = scaled_kinship[scaled_kinship[trait].notna()]
        locations = trait_data['location'].unique()

        expected = 0
        for train_loc in locations:
            target_data = trait_data[trait_data['location'] != train_loc]
            target_lys = target_data.groupby('location_year')['sample.id'].nunique()
            expected += (target_lys >= min_g).sum()

        assert len(results) == expected


# ── run_all_cv_schemes (integration) ────────────────────────────────────────

class TestRunAllCVSchemes:

    def test_all_scheme_keys_present(self, scaled_kinship):
        out = run_all_cv_schemes(scaled_kinship, ['TGW'],
                                 model_class=_QuickRegressor,
                                 k_folds=3, n_iterations=1,
                                 min_genotypes=3)
        assert {'cv1', 'cv2', 'cv0', 'cross_loc', 'all_results'} == set(out.keys())

    def test_all_results_combines_schemes(self, scaled_kinship):
        out = run_all_cv_schemes(scaled_kinship, ['TGW'],
                                 model_class=_QuickRegressor,
                                 k_folds=3, n_iterations=1,
                                 min_genotypes=3)
        all_res = out['all_results']
        if len(all_res) > 0:
            schemes = all_res['cv_scheme'].unique()
            assert len(schemes) >= 2
