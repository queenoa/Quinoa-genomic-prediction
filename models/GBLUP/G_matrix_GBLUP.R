# Load library
library(ASRgenomics)  # Genomic relationship matrix processing

# Load previously created kinship matrix
load("/Users/stansccs/Documents/postdoc/Collaborations/Quinoa-genomic-prediction-code/data/kinship_matrix_VanRaden_auspak_maxmissing20.RData")

# Apply bending to ensure G matrix is positive definite
# eig.tol = 1e-06 is the default tolerance for minimum eigenvalue
G_bending <- G.tuneup(G = kinship_matrix_V2, bend = TRUE)

# Compute inverse of bent G matrix in sparse format for ASReml
# Sparse format reduces memory usage and improves computational speed
Ginv_sparse <- G.inverse(G_bending$Gb, sparse = TRUE)$Ginv

# Save for use by GBLUP.R
save(Ginv_sparse, file = "Ginv_sparse_GBLUP.RData")