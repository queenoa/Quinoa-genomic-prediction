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
    pheno = pd.read_csv('../../data/AUSPAK_phenotypes_GP_input.csv')

    # kinship matrix
    kinship = pd.read_csv('../../data/kinship_matrix_VanRaden_auspak_maxmissing20.csv', index_col=0)

    # Ensure output directory exists
    os.makedirs('model_inputs', exist_ok=True)


    ## merge and prepare model input data

    # Prefix kinship columns with 'K_' so they are distinguishable as features
    kinship.index.name = 'sample.id'
    kinship_df = kinship.reset_index()
    kinship_df.columns = ['sample.id'] + [f'K_{c}' for c in kinship.columns]

    model_input_kinship = pd.merge(kinship_df, pheno, on='sample.id')

    ### One-hot encode
    model_input_final = one_hot_encode(model_input_kinship)

    ### save prepared data
    model_input_final.to_pickle('model_inputs/model_input.pkl')
    print(f"Kinship: {model_input_final.shape[0]} observations, {model_input_final.shape[1]} features")
    print(f"  Saved to: model_inputs/model_input.pkl")
