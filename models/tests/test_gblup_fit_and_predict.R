library(testthat)
library(ASRgenomics)

suppressPackageStartupMessages(
  source(file.path("..", "models", "GBLUP", "GBLUP_utils.R"))
)

# ── Build Ginv_sparse from test marker data ─────────────────────────────────
#
# Mirrors the production pipeline in G_matrix_GBLUP.R:
#   1. Compute VanRaden method 1 G matrix from markers
#   2. Apply bending via ASRgenomics::G.tuneup() for positive definiteness
#   3. Compute sparse inverse via ASRgenomics::G.inverse() for ASReml

MARKER_FILE <- file.path("..", "data", "AUSPAK_test_subset_1k.raw")
PHENO_FILE  <- file.path("..", "data", "AUSPAK_phenotypes_means_BLUEs.csv")

build_test_ginv <- function(marker_file) {
  marker_data <- read.table(marker_file, header = TRUE)
  geno_ids <- as.character(marker_data$IID)
  X <- as.matrix(marker_data[, 7:ncol(marker_data)])
  rownames(X) <- geno_ids

  # Mean-impute missing genotypes
  col_means <- colMeans(X, na.rm = TRUE)
  for (j in seq_len(ncol(X))) {
    na_idx <- is.na(X[, j])
    if (any(na_idx)) X[na_idx, j] <- col_means[j]
  }

  # Remove monomorphic markers
  col_vars <- apply(X, 2, var)
  X <- X[, col_vars > 0]

  # VanRaden method 1
  p <- colMeans(X) / 2
  W <- sweep(X, 2, 2 * p)
  denom <- 2 * sum(p * (1 - p))
  G <- tcrossprod(W) / denom
  rownames(G) <- colnames(G) <- geno_ids

  # Bend for positive definiteness (matches G_matrix_GBLUP.R)
  G_bent <- G.tuneup(G = G, bend = TRUE)$Gb

  # Sparse inverse for ASReml (matches G_matrix_GBLUP.R)
  Ginv_sparse <- G.inverse(G_bent, sparse = TRUE)$Ginv

  Ginv_sparse
}

# ── Cached test data setup ──────────────────────────────────────────────────

test_env <- new.env(parent = emptyenv())

setup_test_data <- function() {
  if (!is.null(test_env$pheno)) return(invisible())

  test_env$Ginv_sparse <- build_test_ginv(MARKER_FILE)

  pheno <- read.csv(PHENO_FILE, stringsAsFactors = FALSE)
  # G.inverse sparse format stores IDs in rowNames attribute, not rownames()
  geno_ids <- attr(test_env$Ginv_sparse, "rowNames")
  pheno <- pheno[pheno$sample.id %in% geno_ids, ]

  # Ensure required columns
  if (!"location_year" %in% names(pheno)) {
    pheno$location_year <- paste(pheno$location, pheno$year, sep = "_")
  }

  test_env$pheno <- pheno
}

# Helper: sample balanced across location-years for multi-environment tests
sample_multi_ly <- function(pheno, n_per_ly = 40) {
  lys <- unique(pheno$location_year)
  do.call(rbind, lapply(lys, function(ly) {
    rows <- pheno[pheno$location_year == ly, ]
    rows[seq_len(min(n_per_ly, nrow(rows))), ]
  }))
}

# ============================================================================
# fit_gblup_and_predict: basic output structure
# ============================================================================

test_that("fit_gblup_and_predict returns dataframe with expected columns", {
  setup_test_data()
  model_data <- sample_multi_ly(test_env$pheno, n_per_ly = 40)
  model_data <- align_genotypes_to_gmatrix(model_data, test_env$Ginv_sparse)

  # Z-score the trait
  model_data <- apply_location_year_scaling(model_data, "PtHt_blue")

  pred <- fit_gblup_and_predict(
    model_data  = model_data,
    trait       = "PtHt_blue",
    Ginv_sparse = test_env$Ginv_sparse
  )

  expect_false(is.null(pred))
  expect_true(is.data.frame(pred))
  expected_cols <- c("sample.id", "location", "year", "location_year",
                     "predicted.value")
  expect_true(all(expected_cols %in% names(pred)),
              info = paste("Missing:", paste(setdiff(expected_cols, names(pred)),
                                              collapse = ", ")))
})

test_that("sample.id, location, year are character (not factor)", {
  setup_test_data()
  model_data <- sample_multi_ly(test_env$pheno, n_per_ly = 40)
  model_data <- align_genotypes_to_gmatrix(model_data, test_env$Ginv_sparse)
  model_data <- apply_location_year_scaling(model_data, "PtHt_blue")

  pred <- fit_gblup_and_predict(model_data, "PtHt_blue", test_env$Ginv_sparse)

  expect_false(is.null(pred))
  expect_true(is.character(pred$sample.id))
  expect_true(is.character(pred$location))
  expect_true(is.character(pred$year))
})

# ============================================================================
# Predictions filter to real location-years
# ============================================================================

test_that("predictions only contain real location-years from input data", {
  setup_test_data()
  model_data <- sample_multi_ly(test_env$pheno, n_per_ly = 40)
  model_data <- align_genotypes_to_gmatrix(model_data, test_env$Ginv_sparse)
  model_data <- apply_location_year_scaling(model_data, "PtHt_blue")

  real_lys <- unique(as.character(model_data$location_year))

  pred <- fit_gblup_and_predict(model_data, "PtHt_blue", test_env$Ginv_sparse)

  expect_false(is.null(pred))
  # All predicted location-years should be in the input data
  expect_true(all(pred$location_year %in% real_lys))
})

# ============================================================================
# Predictions for held-out (NA) rows
# ============================================================================

test_that("predictions exist for held-out genotypes (NA trait values)", {
  setup_test_data()
  model_data <- sample_multi_ly(test_env$pheno, n_per_ly = 40)
  model_data <- align_genotypes_to_gmatrix(model_data, test_env$Ginv_sparse)
  model_data <- apply_location_year_scaling(model_data, "PtHt_blue")

  # Mask some genotypes (CV1-style)
  genos <- unique(as.character(model_data$sample.id))
  set.seed(42)
  test_genos <- sample(genos, min(10, length(genos) %/% 3))
  mask_rows <- model_data$sample.id %in% test_genos
  model_data[["PtHt_blue"]][mask_rows] <- NA

  pred <- fit_gblup_and_predict(model_data, "PtHt_blue", test_env$Ginv_sparse)

  expect_false(is.null(pred))

  # ASReml with na.method(y="include") should predict for ALL genotypes,
  # including those with NA
  pred_genos <- unique(pred$sample.id)
  expect_true(all(test_genos %in% pred_genos),
              info = "Held-out genotypes should have predictions")
})

# ============================================================================
# Error handling: insufficient training data
# ============================================================================

test_that("returns NULL with warning for insufficient training data", {
  setup_test_data()
  model_data <- sample_multi_ly(test_env$pheno, n_per_ly = 10)
  model_data <- align_genotypes_to_gmatrix(model_data, test_env$Ginv_sparse)

  # Set almost all trait values to NA
  n <- nrow(model_data)
  model_data[["PtHt_blue"]][1:(n - 5)] <- NA

  expect_warning(
    pred <- fit_gblup_and_predict(model_data, "PtHt_blue", test_env$Ginv_sparse),
    "Insufficient training data"
  )

  expect_null(pred)
})

# ============================================================================
# Model structure: location fixed + G random + location:year random
# ============================================================================

test_that("predictions have non-trivial variance (model actually fitted)", {
  setup_test_data()
  model_data <- sample_multi_ly(test_env$pheno, n_per_ly = 40)
  model_data <- align_genotypes_to_gmatrix(model_data, test_env$Ginv_sparse)
  model_data <- apply_location_year_scaling(model_data, "PtHt_blue")

  pred <- fit_gblup_and_predict(model_data, "PtHt_blue", test_env$Ginv_sparse)

  expect_false(is.null(pred))
  expect_gt(sd(pred$predicted.value), 0.01,
            label = "Prediction variance should be non-trivial")
})

test_that("predictions span multiple location-years", {
  setup_test_data()
  model_data <- sample_multi_ly(test_env$pheno, n_per_ly = 40)
  model_data <- align_genotypes_to_gmatrix(model_data, test_env$Ginv_sparse)
  model_data <- apply_location_year_scaling(model_data, "PtHt_blue")

  pred <- fit_gblup_and_predict(model_data, "PtHt_blue", test_env$Ginv_sparse)

  expect_false(is.null(pred))
  expect_gt(length(unique(pred$location_year)), 1)
})

# ============================================================================
# location_year is constructed as location + "_" + year
# ============================================================================

test_that("location_year column is paste(location, year, sep='_')", {
  setup_test_data()
  model_data <- sample_multi_ly(test_env$pheno, n_per_ly = 40)
  model_data <- align_genotypes_to_gmatrix(model_data, test_env$Ginv_sparse)
  model_data <- apply_location_year_scaling(model_data, "PtHt_blue")

  pred <- fit_gblup_and_predict(model_data, "PtHt_blue", test_env$Ginv_sparse)

  expect_false(is.null(pred))
  expected_ly <- paste(pred$location, pred$year, sep = "_")
  expect_equal(pred$location_year, expected_ly)
})
