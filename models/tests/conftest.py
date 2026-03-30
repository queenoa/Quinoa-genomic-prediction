"""Shared fixtures for LightGBM pipeline tests.

Uses subsets of real AUSPAK data (created by create_test_fixtures.py).
"""

import sys
import os
import numpy as np
import pandas as pd
import pytest

# Add parent directory to path so we can import the pipeline modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

FIXTURE_DIR = os.path.join(os.path.dirname(__file__), 'fixtures')


@pytest.fixture
def rng():
    """Fixed random state for reproducibility."""
    return np.random.RandomState(42)


@pytest.fixture
def pheno_subset():
    """Real phenotype subset — 100 genotypes, 6 location-years."""
    return pd.read_csv(os.path.join(FIXTURE_DIR, 'pheno_subset.csv'),
                       index_col=0)


@pytest.fixture
def pca_subset():
    """Real PCA subset — 100 genotypes, 551 PCs."""
    return pd.read_csv(os.path.join(FIXTURE_DIR, 'pca_subset.csv'))


@pytest.fixture
def kinship_subset():
    """Real kinship subset — 100 × 100."""
    return pd.read_csv(os.path.join(FIXTURE_DIR, 'kinship_subset.csv'),
                       index_col=0)


@pytest.fixture
def model_input_pc(pheno_subset, pca_subset):
    """Prepared model input with all PCs and one-hot encoding.

    Mirrors the output of Prepare_input_data.py for the all-PC approach.
    """
    from Prepare_input_data import one_hot_encode

    merged = pd.merge(pca_subset, pheno_subset.reset_index(), on='sample.id')
    return one_hot_encode(merged)


@pytest.fixture
def model_input_kinship(pheno_subset, kinship_subset):
    """Prepared model input with kinship features and one-hot encoding."""
    from Prepare_input_data import one_hot_encode

    kinship_subset.index.name = 'sample.id'
    kin_df = kinship_subset.reset_index()
    kin_df.columns = ['sample.id'] + [f'K_{c}' for c in kinship_subset.columns]

    merged = pd.merge(kin_df, pheno_subset.reset_index(), on='sample.id')
    return one_hot_encode(merged)


@pytest.fixture
def predictions_df(rng):
    """Sample predictions DataFrame for evaluate_per_location_year."""
    genos = [f'G{i:03d}' for i in range(1, 21)]  # 20 genotypes
    lys = ['LocA_2020', 'LocA_2021', 'LocB_2020']

    rows = []
    for ly in lys:
        loc = ly.split('_')[0]
        for g in genos:
            obs = rng.randn()
            rows.append({
                'sample.id': g,
                'location_year': ly,
                'location': loc,
                'observed': obs,
                'predicted': obs + rng.randn() * 0.3,
            })

    return pd.DataFrame(rows)
