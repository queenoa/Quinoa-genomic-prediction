library(testthat)
library(ASRgenomics)

# Sourcing GBLUP_per_location_CV1.R is gated by sys.nframe(); inside source()
# from this test (which is itself sourced via test_dir()), sys.nframe() > 0
# so the run-pipeline block at the end is skipped, leaving us with just the
# function definitions (fit_gblup_per_location, fit_single_location_gblup,
# run_per_location_variances, run_cv1_per_location).
old_wd <- setwd(file.path("..", "GBLUP"))
suppressPackageStartupMessages(source("GBLUP_per_location_CV1.R"))
setwd(old_wd)

MARKER_FILE <- file.path("..", "..", "data", "AUSPAK_test_subset_1k.raw")
PHENO_FILE  <- file.path("..", "..", "data", "AUSPAK_phenotypes_GP_input.csv")

# ── Build Ginv_sparse from test marker data ─────────────────────────────────

build_test_ginv <- function(marker_file) {
  marker_data <- read.table(marker_file, header = TRUE)
  geno_ids <- as.character(marker_data$IID)
  X <- as.matrix(marker_data[, 7:ncol(marker_data)])
  rownames(X) <- geno_ids

  col_means <- colMeans(X, na.rm = TRUE)
  for (j in seq_len(ncol(X))) {
    na_idx <- is.na(X[, j])
    if (any(na_idx)) X[na_idx, j] <- col_means[j]
  }

  col_vars <- apply(X, 2, var)
  X <- X[, col_vars > 0]

  p <- colMeans(X) / 2
  W <- sweep(X, 2, 2 * p)
  denom <- 2 * sum(p * (1 - p))
  G <- tcrossprod(W) / denom
  rownames(G) <- colnames(G) <- geno_ids

  G_bent <- G.tuneup(G = G, bend = TRUE)$Gb
  G.inverse(G_bent, sparse = TRUE)$Ginv
}

# ── Cached test data ────────────────────────────────────────────────────────

test_env <- new.env(parent = emptyenv())

setup_test_data <- function() {
  if (!is.null(test_env$pheno)) return(invisible())

  test_env$Ginv_sparse <- build_test_ginv(MARKER_FILE)

  pheno <- read.csv(PHENO_FILE, stringsAsFactors = FALSE)
  geno_ids <- attr(test_env$Ginv_sparse, "rowNames")
  pheno <- pheno[pheno$sample.id %in% geno_ids, ]

  if (!"location_year" %in% names(pheno)) {
    pheno$location_year <- paste(pheno$location, pheno$year, sep = "_")
  }

  test_env$pheno <- pheno
}

TEST_TRAIT <- "PtHt"

# Subset a single location's data (matches what run_cv1_per_location passes
# to fit_gblup_per_location)
prep_loc_data <- function(loc) {
  setup_test_data()
  pheno <- test_env$pheno
  pheno <- pheno[pheno$location == loc, ]
  pheno <- align_genotypes_to_gmatrix(pheno, test_env$Ginv_sparse)
  apply_location_year_scaling(pheno, TEST_TRAIT)
}

# ============================================================================
# fit_gblup_per_location: CV fit on one location's subset
# ============================================================================

test_that("fit_gblup_per_location returns dataframe with expected columns", {
  loc <- "AUS"
  loc_data <- prep_loc_data(loc)
  # Skip if location has too few observations (test subset is small)
  skip_if(sum(!is.na(loc_data[[TEST_TRAIT]])) < 50,
          paste("Not enough data for", loc))

  pred <- fit_gblup_per_location(loc_data, TEST_TRAIT, test_env$Ginv_sparse)

  expect_false(is.null(pred))
  expect_true(is.data.frame(pred))
  expect_true(all(c("sample.id", "location", "year", "location_year",
                     "predicted.value") %in% names(pred)))
})

test_that("fit_gblup_per_location predictions stay within the focal location", {
  loc <- "AUS"
  loc_data <- prep_loc_data(loc)
  skip_if(sum(!is.na(loc_data[[TEST_TRAIT]])) < 50,
          paste("Not enough data for", loc))

  pred <- fit_gblup_per_location(loc_data, TEST_TRAIT, test_env$Ginv_sparse)
  expect_false(is.null(pred))

  expect_true(all(pred$location == loc),
              info = "Predictions should only contain the focal location")

  real_lys <- unique(as.character(loc_data$location_year))
  expect_true(all(pred$location_year %in% real_lys))
})

test_that("fit_gblup_per_location returns NULL with insufficient training data", {
  loc <- "AUS"
  loc_data <- prep_loc_data(loc)

  # Mask almost all observations
  loc_data[[TEST_TRAIT]] <- NA
  loc_data[[TEST_TRAIT]][1:10] <- rnorm(10)

  expect_warning(
    pred <- fit_gblup_per_location(loc_data, TEST_TRAIT, test_env$Ginv_sparse),
    "Insufficient training data"
  )
  expect_null(pred)
})

# ============================================================================
# fit_single_location_gblup: full-data per-location variance components
# ============================================================================

test_that("fit_single_location_gblup returns variance components and h2", {
  loc <- "AUS"
  loc_data <- prep_loc_data(loc)
  skip_if(sum(!is.na(loc_data[[TEST_TRAIT]])) < 30,
          paste("Not enough data for", loc))

  vc <- fit_single_location_gblup(loc_data, TEST_TRAIT, test_env$Ginv_sparse)

  expect_false(is.null(vc))
  expect_true(is.data.frame(vc))
  expect_named(vc, c("n_obs", "genetic_variance", "residual_variance", "h2"))
  expect_gt(vc$n_obs, 0)
  expect_true(vc$residual_variance >= 0)
})

test_that("fit_single_location_gblup h2 is in [0, 1] when both variances positive", {
  loc <- "AUS"
  loc_data <- prep_loc_data(loc)
  skip_if(sum(!is.na(loc_data[[TEST_TRAIT]])) < 30,
          paste("Not enough data for", loc))

  vc <- fit_single_location_gblup(loc_data, TEST_TRAIT, test_env$Ginv_sparse)
  expect_false(is.null(vc))

  if (!is.na(vc$genetic_variance) && vc$genetic_variance >= 0 &&
      vc$residual_variance > 0) {
    expect_gte(vc$h2, 0)
    expect_lte(vc$h2, 1)
  }
})

test_that("fit_single_location_gblup returns NULL with too few observations", {
  loc <- "AUS"
  loc_data <- prep_loc_data(loc)
  loc_data[[TEST_TRAIT]] <- NA
  loc_data[[TEST_TRAIT]][1:5] <- rnorm(5)

  expect_warning(
    vc <- fit_single_location_gblup(loc_data, TEST_TRAIT, test_env$Ginv_sparse),
    "Insufficient data"
  )
  expect_null(vc)
})

# ============================================================================
# run_per_location_variances: one row per (trait, location)
# ============================================================================

test_that("run_per_location_variances returns trait x location rows", {
  setup_test_data()

  vr <- run_per_location_variances(
    pheno_data  = test_env$pheno,
    Ginv_sparse = test_env$Ginv_sparse,
    traits      = TEST_TRAIT,
    apply_zscore = TRUE
  )

  expect_true(is.data.frame(vr))
  expect_true(all(c("trait", "location", "n_obs", "genetic_variance",
                     "residual_variance", "h2") %in% names(vr)))
  # Should have at least one (trait, location) row (assuming enough data)
  expect_gte(nrow(vr), 1)
  expect_true(all(vr$trait == TEST_TRAIT))
})

# ============================================================================
# run_cv1_per_location: full CV1 within each location
# ============================================================================

test_that("run_cv1_per_location returns expected structure", {
  setup_test_data()

  result <- run_cv1_per_location(
    pheno_data    = test_env$pheno,
    Ginv_sparse   = test_env$Ginv_sparse,
    traits        = TEST_TRAIT,
    k_folds       = 2,
    n_iterations  = 1,
    apply_zscore  = TRUE,
    min_genotypes = 5
  )

  expect_true(is.list(result))
  expect_true(all(c("results", "predictions", "summary",
                     "zscore_applied", "cv_scheme") %in% names(result)))
  expect_equal(result$cv_scheme, "CV1_per_location")
  expect_true(result$zscore_applied)
})

test_that("run_cv1_per_location predictions are tagged with train_location", {
  setup_test_data()

  result <- run_cv1_per_location(
    pheno_data    = test_env$pheno,
    Ginv_sparse   = test_env$Ginv_sparse,
    traits        = TEST_TRAIT,
    k_folds       = 2,
    n_iterations  = 1,
    apply_zscore  = TRUE,
    min_genotypes = 5
  )

  expect_true("train_location" %in% names(result$predictions))
  expect_true("train_location" %in% names(result$results))

  if (nrow(result$predictions) > 0) {
    # train_location should match the location column for every row
    # (per-location model predicts only within its own location)
    expect_true(all(result$predictions$train_location ==
                    result$predictions$location))
  }
})

test_that("run_cv1_per_location cv_scheme encodes the location", {
  setup_test_data()

  result <- run_cv1_per_location(
    pheno_data    = test_env$pheno,
    Ginv_sparse   = test_env$Ginv_sparse,
    traits        = TEST_TRAIT,
    k_folds       = 2,
    n_iterations  = 1,
    apply_zscore  = TRUE,
    min_genotypes = 5
  )

  if (nrow(result$results) > 0) {
    expect_true(all(grepl("^CV1_", result$results$cv_scheme)))
  }
})

test_that("run_cv1_per_location uses CV1 seed scheme (1000 + iter)", {
  setup_test_data()

  result <- run_cv1_per_location(
    pheno_data    = test_env$pheno,
    Ginv_sparse   = test_env$Ginv_sparse,
    traits        = TEST_TRAIT,
    k_folds       = 2,
    n_iterations  = 2,
    apply_zscore  = TRUE,
    min_genotypes = 5
  )

  if (nrow(result$results) > 0) {
    for (iter in unique(result$results$iteration)) {
      iter_seeds <- unique(result$results$seed[result$results$iteration == iter])
      expect_equal(iter_seeds, 1000 + iter,
                   info = paste("Iter", iter, "seed should be", 1000 + iter))
    }
  }
})
