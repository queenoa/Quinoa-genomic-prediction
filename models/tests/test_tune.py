"""Tests for tune_LightGBM.py — convert_numpy_types and tune_trait."""

import numpy as np
import pandas as pd
import pytest

from tune_LightGBM import convert_numpy_types, tune_trait, PARAM_DISTRIBUTIONS
from run_LightGBM import get_feature_columns
from LightGBM_utils import apply_location_year_scaling, VALID_TRAITS


# ── convert_numpy_types ─────────────────────────────────────────────────────

class TestConvertNumpyTypes:

    def test_numpy_int(self):
        assert convert_numpy_types(np.int64(5)) == 5
        assert isinstance(convert_numpy_types(np.int64(5)), int)

    def test_numpy_float(self):
        assert convert_numpy_types(np.float64(3.14)) == pytest.approx(3.14)
        assert isinstance(convert_numpy_types(np.float64(3.14)), float)

    def test_numpy_array(self):
        arr = np.array([1, 2, 3])
        assert convert_numpy_types(arr) == [1, 2, 3]

    def test_nested_dict(self):
        d = {'a': np.int64(1), 'b': {'c': np.float64(2.5)}}
        result = convert_numpy_types(d)
        assert result == {'a': 1, 'b': {'c': 2.5}}
        assert isinstance(result['a'], int)
        assert isinstance(result['b']['c'], float)

    def test_plain_python_types_unchanged(self):
        assert convert_numpy_types(42) == 42
        assert convert_numpy_types('hello') == 'hello'
        assert convert_numpy_types(3.14) == 3.14


# ── tune_trait ──────────────────────────────────────────────────────────────

class TestTuneTrait:

    @pytest.fixture
    def tuning_data(self, model_input_kinship):
        """Scaled real data with feature columns for tuning."""
        traits = [t for t in VALID_TRAITS if t in model_input_kinship.columns]
        scaled = apply_location_year_scaling(model_input_kinship, traits)
        feature_cols = get_feature_columns(scaled)
        return scaled, feature_cols

    def test_returns_best_params(self, tuning_data):
        df, feature_cols = tuning_data
        best, score = tune_trait(df, 'TGW', feature_cols)
        assert best is not None
        assert isinstance(best, dict)
        assert any(k in best for k in PARAM_DISTRIBUTIONS)

    def test_returns_negative_mse_score(self, tuning_data):
        df, feature_cols = tuning_data
        best, score = tune_trait(df, 'TGW', feature_cols)
        assert score < 0  # neg_mean_squared_error is always negative

    def test_skips_trait_with_few_observations(self, tuning_data):
        df, feature_cols = tuning_data
        # Keep only 50 rows (< 100 threshold)
        small_df = df.head(50).copy()
        best, score = tune_trait(small_df, 'TGW', feature_cols)
        assert best is None
        assert score is None
