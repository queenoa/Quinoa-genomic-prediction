"""Tests for LightGBM_utils.py — constants, metrics, scaling, evaluation."""

import numpy as np
import pandas as pd
import pytest

from LightGBM_utils import (
    VALID_TRAITS, LOWER_IS_BETTER_TRAITS, DEFAULT_PARAMS,
    fetch_model_params, calculate_ndcg, evaluate_predictions,
    apply_location_year_scaling, evaluate_per_location_year,
    summarise_cv_results,
)


# ── Constants ───────────────────────────────────────────────────────────────

class TestConstants:

    def test_valid_traits_count(self):
        assert len(VALID_TRAITS) == 7

    def test_lower_is_better_subset(self):
        assert set(LOWER_IS_BETTER_TRAITS).issubset(set(VALID_TRAITS))

    def test_lower_is_better_contains_expected(self):
        assert 'DTF_blue' in LOWER_IS_BETTER_TRAITS
        assert 'DTH_blue' in LOWER_IS_BETTER_TRAITS
        assert 'PtHt_blue' in LOWER_IS_BETTER_TRAITS

    def test_default_params_keys(self):
        assert set(DEFAULT_PARAMS.keys()) == {'max_depth', 'learning_rate', 'n_estimators'}


# ── fetch_model_params ──────────────────────────────────────────────────────

class TestFetchModelParams:

    def test_none_returns_defaults(self):
        assert fetch_model_params(None, 'DTF_blue') == DEFAULT_PARAMS

    def test_flat_dict_returned_as_is(self):
        flat = {'max_depth': 5, 'learning_rate': 0.1, 'n_estimators': 200}
        assert fetch_model_params(flat, 'DTF_blue') == flat

    def test_per_trait_dict(self):
        per_trait = {
            'DTF_blue': {'max_depth': 4, 'learning_rate': 0.03, 'n_estimators': 1000},
            'TGW_blue': {'max_depth': 6, 'learning_rate': 0.1, 'n_estimators': 200},
        }
        result = fetch_model_params(per_trait, 'DTF_blue')
        assert result == per_trait['DTF_blue']

    def test_per_trait_dict_missing_trait_returns_whole_dict(self):
        per_trait = {
            'DTF_blue': {'max_depth': 4},
        }
        # SdLen_blue not in dict → returns the outer dict itself (flat fallback)
        result = fetch_model_params(per_trait, 'SdLen_blue')
        assert result == per_trait


# ── calculate_ndcg ──────────────────────────────────────────────────────────

class TestCalculateNDCG:

    def test_perfect_ranking(self):
        y = np.array([5.0, 4.0, 3.0, 2.0, 1.0])
        score = calculate_ndcg(y, y, k=5)
        assert score == pytest.approx(1.0)

    def test_returns_between_0_and_1(self):
        rng = np.random.RandomState(7)
        y_true = rng.randn(20)
        y_pred = rng.randn(20)
        score = calculate_ndcg(y_true, y_pred, k=10)
        assert 0.0 <= score <= 1.0

    def test_lower_is_better_flips_sign(self):
        y_true = np.array([1.0, 2.0, 3.0, 4.0, 5.0])  # lower=better → best is 1
        y_pred = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        score = calculate_ndcg(y_true, y_pred, k=5, lower_is_better=True)
        assert score == pytest.approx(1.0)

    def test_k_clamped_to_array_length(self):
        y = np.array([3.0, 1.0, 2.0])
        score = calculate_ndcg(y, y, k=100)
        assert score == pytest.approx(1.0)

    def test_negative_values_handled(self):
        y_true = np.array([-3.0, -1.0, 0.0, 2.0])
        y_pred = np.array([-2.5, -0.5, 0.5, 1.5])
        score = calculate_ndcg(y_true, y_pred, k=4)
        assert 0.0 <= score <= 1.0


# ── evaluate_predictions ────────────────────────────────────────────────────

class TestEvaluatePredictions:

    def test_perfect_prediction(self):
        y = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
        result = evaluate_predictions(y, y, trait_name='TGW_blue')
        assert result['pearson'] == pytest.approx(1.0)
        assert result['spearman'] == pytest.approx(1.0)
        assert result['ndcg_at_10'] == pytest.approx(1.0)

    def test_returns_nan_for_single_value(self):
        result = evaluate_predictions(np.array([1.0]), np.array([2.0]))
        assert np.isnan(result['pearson'])
        assert np.isnan(result['spearman'])
        assert np.isnan(result['ndcg_at_10'])

    def test_output_keys(self):
        y = np.arange(10, dtype=float)
        result = evaluate_predictions(y, y + 0.1)
        assert set(result.keys()) == {'pearson', 'spearman', 'ndcg_at_10'}

    def test_lower_is_better_trait_recognised(self):
        y_true = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
        result = evaluate_predictions(y_true, y_true, trait_name='DTF_blue')
        # Perfect ranking should still give 1.0 regardless of direction
        assert result['ndcg_at_10'] == pytest.approx(1.0)

    def test_no_trait_name_defaults_higher_is_better(self):
        y = np.arange(10, dtype=float)
        result = evaluate_predictions(y, y, trait_name=None)
        assert result['ndcg_at_10'] == pytest.approx(1.0)


# ── apply_location_year_scaling ─────────────────────────────────────────────

class TestApplyLocationYearScaling:

    def test_output_has_zero_mean_unit_std(self, model_input_pc):
        traits = ['DTF_blue', 'TGW_blue']
        scaled = apply_location_year_scaling(model_input_pc, traits)

        for trait in traits:
            for ly in scaled['location_year'].unique():
                vals = scaled.loc[
                    (scaled['location_year'] == ly) & scaled[trait].notna(), trait
                ]
                if len(vals) > 1:
                    assert vals.mean() == pytest.approx(0.0, abs=1e-10)
                    assert vals.std() == pytest.approx(1.0, abs=1e-10)

    def test_does_not_modify_original(self, model_input_pc):
        original = model_input_pc.copy()
        apply_location_year_scaling(model_input_pc, ['TGW_blue'])
        pd.testing.assert_frame_equal(model_input_pc, original)

    def test_skips_trait_with_few_observations(self):
        """Traits with <50 non-null values should be left unchanged."""
        df = pd.DataFrame({
            'sample.id': [f'G{i}' for i in range(10)],
            'location_year': ['L1'] * 10,
            'small_trait': np.arange(10, dtype=float),
        })
        scaled = apply_location_year_scaling(df, ['small_trait'])
        np.testing.assert_array_equal(
            scaled['small_trait'].values, df['small_trait'].values
        )

    def test_zero_std_location_year_uses_1(self):
        """If all values in a location-year are equal, std=0 → use 1."""
        n = 60
        df = pd.DataFrame({
            'sample.id': [f'G{i}' for i in range(n)],
            'location_year': ['LY1'] * n,
            'trait': [5.0] * n,
        })
        scaled = apply_location_year_scaling(df, ['trait'])
        # (5 - 5) / 1 = 0
        assert (scaled['trait'] == 0.0).all()

    def test_nans_preserved(self, model_input_pc):
        nan_before = model_input_pc['DTF_blue'].isna().sum()
        scaled = apply_location_year_scaling(model_input_pc, ['DTF_blue'])
        nan_after = scaled['DTF_blue'].isna().sum()
        assert nan_after == nan_before


# ── evaluate_per_location_year ──────────────────────────────────────────────

class TestEvaluatePerLocationYear:

    def test_returns_one_row_per_location_year(self, predictions_df):
        result = evaluate_per_location_year(predictions_df, 'TGW_blue',
                                            min_genotypes=10)
        assert len(result) == predictions_df['location_year'].nunique()

    def test_skips_location_year_below_min_genotypes(self, predictions_df):
        # Keep only 5 genotypes in one location-year
        small = predictions_df[
            (predictions_df['location_year'] == 'LocA_2020') &
            (predictions_df['sample.id'].isin(['G001', 'G002', 'G003', 'G004', 'G005']))
        ]
        rest = predictions_df[predictions_df['location_year'] != 'LocA_2020']
        df = pd.concat([small, rest], ignore_index=True)

        result = evaluate_per_location_year(df, 'TGW_blue', min_genotypes=10)
        assert 'LocA_2020' not in result['location_year'].values

    def test_output_columns(self, predictions_df):
        result = evaluate_per_location_year(predictions_df, 'TGW_blue')
        expected_cols = {'trait', 'location_year', 'location', 'pearson',
                         'spearman', 'ndcg_at_10', 'n_test_genotypes'}
        assert expected_cols == set(result.columns)

    def test_empty_when_all_below_threshold(self):
        """All location-years have fewer genotypes than min → empty result."""
        df = pd.DataFrame({
            'sample.id': ['G1', 'G2'],
            'location_year': ['LY1', 'LY1'],
            'location': ['L1', 'L1'],
            'observed': [1.0, 2.0],
            'predicted': [1.1, 2.1],
        })
        result = evaluate_per_location_year(df, 'TGW_blue', min_genotypes=10)
        assert len(result) == 0


# ── summarise_cv_results ────────────────────────────────────────────────────

class TestSummariseCVResults:

    def test_summary_has_expected_columns(self, predictions_df):
        result = evaluate_per_location_year(predictions_df, 'TGW_blue')
        summary = summarise_cv_results(result, 'TestScheme')
        expected = {'trait', 'location', 'pearson_mean', 'pearson_std',
                    'pearson_min', 'pearson_max', 'spearman_mean', 'spearman_std',
                    'spearman_min', 'spearman_max', 'ndcg_at_10_mean',
                    'ndcg_at_10_std', 'ndcg_at_10_min', 'ndcg_at_10_max',
                    'n_evaluations', 'mean_n_test_genotypes'}
        assert expected.issubset(set(summary.columns))

    def test_empty_input_returns_empty(self):
        summary = summarise_cv_results(pd.DataFrame(), 'Empty')
        assert len(summary) == 0
