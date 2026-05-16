"""Tests for Prepare_input_data.py — one_hot_encode function.

Uses real data subsets loaded via conftest fixtures. Production pipeline is
kinship-only (commit c97e97a dropped PC features).
"""

import numpy as np
import pandas as pd
import pytest

from Prepare_input_data import one_hot_encode


class TestOneHotEncode:

    @pytest.fixture
    def raw_merge(self, pheno_subset, kinship_subset):
        """Merged kinship + pheno before one-hot encoding (mirrors Prepare_input_data.py)."""
        kinship_subset.index.name = 'sample.id'
        kin_df = kinship_subset.reset_index()
        kin_df.columns = ['sample.id'] + [f'K_{c}' for c in kinship_subset.columns]
        return pd.merge(kin_df, pheno_subset, on='sample.id')

    def test_year_column_dropped(self, raw_merge):
        result = one_hot_encode(raw_merge)
        assert 'year' not in result.columns

    def test_one_hot_location_columns_created(self, raw_merge):
        result = one_hot_encode(raw_merge)
        loc_cols = [c for c in result.columns if c.startswith('location_')
                    and not c.startswith('location_year')]
        n_locations = raw_merge['location'].nunique()
        assert len(loc_cols) == n_locations

    def test_one_hot_location_year_columns_created(self, raw_merge):
        result = one_hot_encode(raw_merge)
        ly_cols = [c for c in result.columns if c.startswith('location_year_')]
        n_ly = raw_merge['location_year'].nunique()
        assert len(ly_cols) == n_ly

    def test_one_hot_values_are_binary(self, raw_merge):
        result = one_hot_encode(raw_merge)
        oh_cols = [c for c in result.columns if c.startswith('location_')
                   and c not in ('location', 'location_year')]
        for col in oh_cols:
            assert set(result[col].unique()).issubset({0.0, 1.0})

    def test_original_columns_preserved(self, raw_merge):
        result = one_hot_encode(raw_merge)
        assert 'sample.id' in result.columns
        # First kinship column should still be present after one-hot encoding
        assert any(c.startswith('K_') for c in result.columns)
        assert 'DTF' in result.columns
        assert 'location' in result.columns
        assert 'location_year' in result.columns

    def test_row_count_unchanged(self, raw_merge):
        result = one_hot_encode(raw_merge)
        assert len(result) == len(raw_merge)

    def test_one_hot_location_sums_to_one(self, raw_merge):
        """Each row should have exactly one 1.0 across location one-hot cols."""
        result = one_hot_encode(raw_merge)
        loc_cols = [c for c in result.columns if c.startswith('location_')
                    and not c.startswith('location_year')]
        row_sums = result[loc_cols].sum(axis=1)
        assert (row_sums == 1.0).all()

    def test_one_hot_location_year_sums_to_one(self, raw_merge):
        result = one_hot_encode(raw_merge)
        ly_cols = [c for c in result.columns if c.startswith('location_year_')]
        row_sums = result[ly_cols].sum(axis=1)
        assert (row_sums == 1.0).all()
