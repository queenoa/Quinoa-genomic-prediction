# Load libraries
import os
import pandas as pd
from sklearn.preprocessing import OneHotEncoder


### ONE-HOT ENCODE CATEGORICAL VARIABLES

def one_hot_encode(model_input):
    """One-hot encode location and location_year, drop year column."""
    cat_encoder = OneHotEncoder(sparse_output=False)
    encoded_features = cat_encoder.fit_transform(model_input[['location', 'location_year']])
    encoded_df = pd.DataFrame(
        encoded_features,
        columns=cat_encoder.get_feature_names_out()
    )

    model_input_final = pd.concat([model_input.reset_index(drop=True), encoded_df], axis=1)
    model_input_final = model_input_final.drop(columns=['year'])
    return model_input_final


if __name__ == '__main__':

    ### load input data

    # phenotype data
    pheno = pd.read_csv('../../data/AUSPAK_phenotypes_means_BLUEs.csv')

    # Drop means columns — only BLUEs are used as traits
    mean_cols = [c for c in pheno.columns if c.endswith('_mean')]
    pheno = pheno.drop(columns=mean_cols)

    # load PCA data
    pca = pd.read_csv('../../data/AUSPAK_PCs_all.csv')

    # kinship matrix
    kinship = pd.read_csv('../../data/kinship_matrix_VanRaden_auspak_maxmissing20.csv', index_col=0)

    # Ensure output directory exists
    os.makedirs('model_inputs', exist_ok=True)


    ## merge and prepare model input data

    # ── Approach 1: All PCs (551) ────────────────────────────────────────────
    model_input_allpc = pd.merge(pca, pheno, on='sample.id')

    # ── Approach 2: 25 PCs ──────────────────────────────────────────────────
    pc_25_cols = ['sample.id'] + [f'PC{i}' for i in range(1, 26)]
    pca_25 = pca[pc_25_cols]
    model_input_25pc = pd.merge(pca_25, pheno, on='sample.id')

    # ── Approach 3: Kinship matrix ──────────────────────────────────────────
    # Prefix kinship columns with 'K_' so they are distinguishable as features
    kinship.index.name = 'sample.id'
    kinship_df = kinship.reset_index()
    kinship_df.columns = ['sample.id'] + [f'K_{c}' for c in kinship.columns]

    model_input_kinship = pd.merge(kinship_df, pheno, on='sample.id')

    ### One-hot encode all three approaches
    model_input_allpc_final = one_hot_encode(model_input_allpc)
    model_input_25pc_final = one_hot_encode(model_input_25pc)
    model_input_kinship_final = one_hot_encode(model_input_kinship)

    ### save prepared data
    model_input_allpc_final.to_pickle('model_inputs/model_input.pkl')
    print(f"All PCs: {model_input_allpc_final.shape[0]} observations, {model_input_allpc_final.shape[1]} features")
    print(f"  Saved to: model_inputs/model_input.pkl")

    model_input_25pc_final.to_pickle('model_inputs/model_input_25pc.pkl')
    print(f"25 PCs: {model_input_25pc_final.shape[0]} observations, {model_input_25pc_final.shape[1]} features")
    print(f"  Saved to: model_inputs/model_input_25pc.pkl")

    model_input_kinship_final.to_pickle('model_inputs/model_input_kinship.pkl')
    print(f"Kinship: {model_input_kinship_final.shape[0]} observations, {model_input_kinship_final.shape[1]} features")
    print(f"  Saved to: model_inputs/model_input_kinship.pkl")
