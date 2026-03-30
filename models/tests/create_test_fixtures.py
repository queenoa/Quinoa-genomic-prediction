"""Create small test fixture files from the real AUSPAK data.

Run once to produce CSV/pickle subsets used by the test suite.
Takes ~20 genotypes so tests run fast but use real data structure.

Usage:
    python tests/create_test_fixtures.py
"""

import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', '..', 'data')
FIXTURE_DIR = os.path.join(os.path.dirname(__file__), 'fixtures')
os.makedirs(FIXTURE_DIR, exist_ok=True)

N_GENOTYPES = 100
SEED = 42

# ── Load full data ──────────────────────────────────────────────────────────

pheno = pd.read_csv(os.path.join(DATA_DIR, 'AUSPAK_phenotypes_means_BLUEs.csv'),
                     index_col=0)
pca = pd.read_csv(os.path.join(DATA_DIR, 'AUSPAK_PCs_all.csv'))
kinship = pd.read_csv(
    os.path.join(DATA_DIR, 'kinship_matrix_VanRaden_auspak_maxmissing20.csv'),
    index_col=0,
)

# ── Sample genotypes ────────────────────────────────────────────────────────

rng = np.random.RandomState(SEED)
all_genos = pheno['sample.id'].unique()
subset_genos = rng.choice(all_genos, size=N_GENOTYPES, replace=False)

# ── Subset phenotype ────────────────────────────────────────────────────────

pheno_sub = pheno[pheno['sample.id'].isin(subset_genos)].copy()
pheno_sub.to_csv(os.path.join(FIXTURE_DIR, 'pheno_subset.csv'))
print(f"Phenotype subset: {pheno_sub.shape[0]} rows, "
      f"{pheno_sub['sample.id'].nunique()} genotypes, "
      f"{pheno_sub.index.nunique()} location-years")

# ── Subset PCA ──────────────────────────────────────────────────────────────

pca_sub = pca[pca['sample.id'].isin(subset_genos)].copy()
pca_sub.to_csv(os.path.join(FIXTURE_DIR, 'pca_subset.csv'), index=False)
print(f"PCA subset: {pca_sub.shape}")

# ── Subset kinship ──────────────────────────────────────────────────────────

kin_sub = kinship.loc[
    kinship.index.isin(subset_genos),
    kinship.columns.isin(subset_genos)
].copy()
kin_sub.to_csv(os.path.join(FIXTURE_DIR, 'kinship_subset.csv'))
print(f"Kinship subset: {kin_sub.shape}")

print(f"\nFixtures saved to {FIXTURE_DIR}/")
