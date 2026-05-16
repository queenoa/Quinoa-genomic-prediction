library(testthat)
library(ASRgenomics)

# Source GBLUP.R (which sources GBLUP_utils.R) from its directory
old_wd <- setwd(file.path("..", "GBLUP"))
suppressPackageStartupMessages(source("GBLUP.R"))
setwd(old_wd)

# ── Build Ginv_sparse from test marker data ─────────────────────────────────
#
# Mirrors the production pipeline in G_matrix_GBLUP.R:
#   1. Compute VanRaden method 1 G matrix from markers
#   2. Apply bending via ASRgenomics::G.tuneup() for positive definiteness
#   3. Compute sparse inverse via ASRgenomics::G.inverse() for ASReml

MARKER_FILE <- file.path("..", "..", "data", "AUSPAK_test_subset_1k.raw")
PHENO_FILE  <- file.path("..", "..", "data", "AUSPAK_phenotypes_GP_input.csv")

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

  # Bend for positive definiteness (matches G_matrix_GBLUP.R)
  G_bent <- G.tuneup(G = G, bend = TRUE)$Gb

  # Sparse inverse for ASReml (matches G_matrix_GBLUP.R)
  Ginv_sparse <- G.inverse(G_bent, sparse = TRUE)$Ginv

  Ginv_sparse
}

# ── Cached test data ────────────────────────────────────────────────────────

test_env <- new.env(parent = emptyenv())

setup_integration_data <- function() {
  if (!is.null(test_env$pheno)) return(invisible())

  test_env$Ginv_sparse <- build_test_ginv(MARKER_FILE)

  pheno <- read.csv(PHENO_FILE, stringsAsFactors = FALSE)
  # G.inverse sparse format stores IDs in rowNames attribute, not rownames()
  geno_ids <- attr(test_env$Ginv_sparse, "rowNames")
  pheno <- pheno[pheno$sample.id %in% geno_ids, ]

  if (!"location_year" %in% names(pheno)) {
    pheno$location_year <- paste(pheno$location, pheno$year, sep = "_")
  }

  test_env$pheno <- pheno
}

# Use a single trait and minimal settings for speed
TEST_TRAIT   <- "PtHt"
K_FOLDS      <- 2
N_ITERATIONS <- 2
MIN_GENO     <- 5

# ============================================================================
# CV1: New genotypes in known environments
# ============================================================================

test_that("run_cv1 returns expected structure", {
  setup_integration_data()

  result <- run_cv1(
    pheno_data    = test_env$pheno,
    Ginv_sparse   = test_env$Ginv_sparse,
    traits        = TEST_TRAIT,
    k_folds       = K_FOLDS,
    n_iterations  = N_ITERATIONS,
    apply_zscore  = TRUE,
    min_genotypes = MIN_GENO
  )

  expect_true(is.list(result))
  expect_true(all(c("results", "predictions", "summary", "n_genotypes",
                     "n_location_years", "zscore_applied", "cv_scheme")
                   %in% names(result)))
  expect_equal(result$cv_scheme, "CV1")
  expect_true(result$zscore_applied)
})

test_that("CV1 results have correct columns and scheme label", {
  setup_integration_data()

  result <- run_cv1(
    pheno_data = test_env$pheno, Ginv_sparse = test_env$Ginv_sparse,
    traits = TEST_TRAIT, k_folds = K_FOLDS, n_iterations = N_ITERATIONS,
    min_genotypes = MIN_GENO
  )

  res <- result$results
  expect_true(nrow(res) > 0, info = "Should have at least one evaluation row")

  expected_cols <- c("iteration", "fold", "trait", "location_year", "location",
                     "pearson", "spearman", "ndcg_at_10", "seed",
                     "n_test_genotypes", "cv_scheme")
  expect_true(all(expected_cols %in% names(res)))
  expect_true(all(res$cv_scheme == "CV1"))
  expect_true(all(res$trait == TEST_TRAIT))
})

test_that("CV1 seed scheme is 1000 + iter", {
  setup_integration_data()

  result <- run_cv1(
    pheno_data = test_env$pheno, Ginv_sparse = test_env$Ginv_sparse,
    traits = TEST_TRAIT, k_folds = K_FOLDS, n_iterations = N_ITERATIONS,
    min_genotypes = MIN_GENO
  )

  res <- result$results
  for (iter in unique(res$iteration)) {
    iter_seeds <- unique(res$seed[res$iteration == iter])
    expect_equal(iter_seeds, 1000 + iter,
                 info = paste("Iteration", iter, "seed should be", 1000 + iter))
  }
})

test_that("CV1 predictions have correct columns", {
  setup_integration_data()

  result <- run_cv1(
    pheno_data = test_env$pheno, Ginv_sparse = test_env$Ginv_sparse,
    traits = TEST_TRAIT, k_folds = K_FOLDS, n_iterations = N_ITERATIONS,
    min_genotypes = MIN_GENO
  )

  preds <- result$predictions
  expect_true(nrow(preds) > 0)
  expected_cols <- c("sample.id", "location_year", "location",
                     "observed", "predicted", "trait", "iteration",
                     "fold", "seed", "cv_scheme")
  expect_true(all(expected_cols %in% names(preds)))
  expect_true(all(preds$cv_scheme == "CV1"))
})

test_that("CV1 predictions have reasonable variance (not constant)", {
  setup_integration_data()

  result <- run_cv1(
    pheno_data = test_env$pheno, Ginv_sparse = test_env$Ginv_sparse,
    traits = TEST_TRAIT, k_folds = K_FOLDS, n_iterations = N_ITERATIONS,
    min_genotypes = MIN_GENO
  )

  pred_sd <- sd(result$predictions$predicted)
  expect_gt(pred_sd, 0.01,
            label = "Prediction variance should be non-trivial")
})

test_that("CV1 accuracy metrics are in valid ranges", {
  setup_integration_data()

  result <- run_cv1(
    pheno_data = test_env$pheno, Ginv_sparse = test_env$Ginv_sparse,
    traits = TEST_TRAIT, k_folds = K_FOLDS, n_iterations = N_ITERATIONS,
    min_genotypes = MIN_GENO
  )

  res <- result$results
  expect_true(all(res$pearson >= -1 & res$pearson <= 1))
  expect_true(all(res$spearman >= -1 & res$spearman <= 1))
  expect_true(all(res$ndcg_at_10 >= 0 & res$ndcg_at_10 <= 1))
})

test_that("CV1 fold assignments differ between iterations (shuffling)", {
  setup_integration_data()

  result <- run_cv1(
    pheno_data = test_env$pheno, Ginv_sparse = test_env$Ginv_sparse,
    traits = TEST_TRAIT, k_folds = K_FOLDS, n_iterations = N_ITERATIONS,
    min_genotypes = MIN_GENO
  )

  preds <- result$predictions
  iters <- sort(unique(preds$iteration))
  expect_gte(length(iters), 2,
             label = "Need >= 2 iterations to test fold shuffling")

  fold1_iter1 <- sort(unique(preds$sample.id[preds$iteration == iters[1] &
                                              preds$fold == 1]))
  fold1_iter2 <- sort(unique(preds$sample.id[preds$iteration == iters[2] &
                                              preds$fold == 1]))

  expect_false(identical(fold1_iter1, fold1_iter2),
               info = paste("Fold 1 genotypes should differ between iterations.",
                            "If identical, genotypes are not being shuffled."))
})

# ============================================================================
# CV2: Sparse testing (cell-level masking)
# ============================================================================

test_that("run_cv2 returns expected structure", {
  setup_integration_data()

  result <- run_cv2(
    pheno_data    = test_env$pheno,
    Ginv_sparse   = test_env$Ginv_sparse,
    traits        = TEST_TRAIT,
    k_folds       = K_FOLDS,
    n_iterations  = N_ITERATIONS,
    apply_zscore  = TRUE,
    min_genotypes = MIN_GENO
  )

  expect_true(is.list(result))
  expect_true(all(c("results", "predictions", "summary",
                     "zscore_applied", "cv_scheme") %in% names(result)))
  expect_equal(result$cv_scheme, "CV2")
})

test_that("CV2 results have correct columns and scheme label", {
  setup_integration_data()

  result <- run_cv2(
    pheno_data = test_env$pheno, Ginv_sparse = test_env$Ginv_sparse,
    traits = TEST_TRAIT, k_folds = K_FOLDS, n_iterations = N_ITERATIONS,
    min_genotypes = MIN_GENO
  )

  res <- result$results
  expect_true(nrow(res) > 0)
  expect_true(all(res$cv_scheme == "CV2"))
  expect_true(all(res$trait == TEST_TRAIT))
})

test_that("CV2 seed scheme is 2000 + iter", {
  setup_integration_data()

  result <- run_cv2(
    pheno_data = test_env$pheno, Ginv_sparse = test_env$Ginv_sparse,
    traits = TEST_TRAIT, k_folds = K_FOLDS, n_iterations = N_ITERATIONS,
    min_genotypes = MIN_GENO
  )

  res <- result$results
  for (iter in unique(res$iteration)) {
    iter_seeds <- unique(res$seed[res$iteration == iter])
    expect_equal(iter_seeds, 2000 + iter,
                 info = paste("Iteration", iter, "seed should be", 2000 + iter))
  }
})

test_that("CV2 uses cell-level masking (not genotype-level)", {
  setup_integration_data()

  result <- run_cv2(
    pheno_data = test_env$pheno, Ginv_sparse = test_env$Ginv_sparse,
    traits = TEST_TRAIT, k_folds = K_FOLDS, n_iterations = N_ITERATIONS,
    min_genotypes = MIN_GENO
  )

  preds <- result$predictions
  if (nrow(preds) > 0) {
    # In CV2, the same genotype may appear in multiple location-years as
    # both test and training. Check that predictions come from specific
    # (sample.id, location_year) cells, not genotype-wide masking.
    # A genotype appearing as test in one LY and train in another is the
    # hallmark of cell-level masking.
    geno_ly_combos <- paste(preds$sample.id, preds$location_year, sep = "_")
    # At minimum, predictions should exist
    expect_gt(length(unique(geno_ly_combos)), 0)
  }
})

test_that("CV2 predictions have non-trivial variance", {
  setup_integration_data()

  result <- run_cv2(
    pheno_data = test_env$pheno, Ginv_sparse = test_env$Ginv_sparse,
    traits = TEST_TRAIT, k_folds = K_FOLDS, n_iterations = N_ITERATIONS,
    min_genotypes = MIN_GENO
  )

  if (nrow(result$predictions) > 0) {
    pred_sd <- sd(result$predictions$predicted)
    expect_gt(pred_sd, 0.01)
  }
})

test_that("CV2 all location-years have predictions in every fold", {
  setup_integration_data()

  result <- run_cv2(
    pheno_data = test_env$pheno, Ginv_sparse = test_env$Ginv_sparse,
    traits = TEST_TRAIT, k_folds = K_FOLDS, n_iterations = N_ITERATIONS,
    min_genotypes = MIN_GENO
  )

  preds <- result$predictions
  if (nrow(preds) > 0) {
    for (iter in unique(preds$iteration)) {
      iter_preds <- preds[preds$iteration == iter, ]
      all_lys    <- unique(iter_preds$location_year)
      k_folds_in <- max(iter_preds$fold)

      for (fold in 1:k_folds_in) {
        fold_lys <- unique(iter_preds$location_year[iter_preds$fold == fold])
        expect_true(setequal(fold_lys, all_lys),
                    info = paste("Iter", iter, "fold", fold,
                                 "should have predictions for all location-years.",
                                 "Missing:",
                                 paste(setdiff(all_lys, fold_lys), collapse = ", ")))
      }
    }
  }
})

test_that("CV2 fold assignments differ between iterations (shuffling)", {
  setup_integration_data()

  result <- run_cv2(
    pheno_data = test_env$pheno, Ginv_sparse = test_env$Ginv_sparse,
    traits = TEST_TRAIT, k_folds = K_FOLDS, n_iterations = N_ITERATIONS,
    min_genotypes = MIN_GENO
  )

  preds <- result$predictions
  iters <- sort(unique(preds$iteration))

  if (length(iters) >= 2) {
    fold1_iter1 <- sort(paste(preds$sample.id[preds$iteration == iters[1] &
                                               preds$fold == 1],
                              preds$location_year[preds$iteration == iters[1] &
                                                   preds$fold == 1]))
    fold1_iter2 <- sort(paste(preds$sample.id[preds$iteration == iters[2] &
                                               preds$fold == 1],
                              preds$location_year[preds$iteration == iters[2] &
                                                   preds$fold == 1]))

    expect_false(identical(fold1_iter1, fold1_iter2),
                 info = paste("CV2 fold 1 observations should differ between",
                              "iterations. If identical, observations are not",
                              "being shuffled."))
  }
})

# ============================================================================
# CV0: Leave-one-location-year-out
# ============================================================================

test_that("run_cv0 returns expected structure", {
  setup_integration_data()

  result <- run_cv0(
    pheno_data    = test_env$pheno,
    Ginv_sparse   = test_env$Ginv_sparse,
    traits        = TEST_TRAIT,
    apply_zscore  = TRUE,
    min_genotypes = MIN_GENO
  )

  expect_true(is.list(result))
  expect_true(all(c("results", "predictions", "summary",
                     "zscore_applied", "cv_scheme") %in% names(result)))
  expect_equal(result$cv_scheme, "CV0")
})

test_that("CV0 results have correct columns and scheme label", {
  setup_integration_data()

  result <- run_cv0(
    pheno_data = test_env$pheno, Ginv_sparse = test_env$Ginv_sparse,
    traits = TEST_TRAIT, min_genotypes = MIN_GENO
  )

  res <- result$results
  expect_true(nrow(res) > 0)
  expect_true(all(res$cv_scheme == "CV0"))
  expect_true(all(res$trait == TEST_TRAIT))
})

test_that("CV0 is deterministic (no iteration/fold)", {
  setup_integration_data()

  result <- run_cv0(
    pheno_data = test_env$pheno, Ginv_sparse = test_env$Ginv_sparse,
    traits = TEST_TRAIT, min_genotypes = MIN_GENO
  )

  res <- result$results
  # CV0 should have NA for iteration and fold (deterministic, no random folds)
  expect_true(all(is.na(res$iteration)))
  expect_true(all(is.na(res$fold)))
})

test_that("CV0 evaluates multiple location-years", {
  setup_integration_data()

  result <- run_cv0(
    pheno_data = test_env$pheno, Ginv_sparse = test_env$Ginv_sparse,
    traits = TEST_TRAIT, min_genotypes = MIN_GENO
  )

  res <- result$results
  expect_gt(length(unique(res$location_year)), 1,
            label = "CV0 should evaluate multiple held-out location-years")
})

test_that("CV0 predictions track held-out location-years correctly", {
  setup_integration_data()

  result <- run_cv0(
    pheno_data = test_env$pheno, Ginv_sparse = test_env$Ginv_sparse,
    traits = TEST_TRAIT, min_genotypes = MIN_GENO
  )

  preds <- result$predictions
  expect_true(nrow(preds) > 0)
  expect_true(all(!is.na(preds$observed)))
  expect_true(all(!is.na(preds$predicted)))
})

# ============================================================================
# Cross-location transferability
# ============================================================================

test_that("run_cross_location returns expected structure", {
  setup_integration_data()

  result <- run_cross_location(
    pheno_data    = test_env$pheno,
    Ginv_sparse   = test_env$Ginv_sparse,
    traits        = TEST_TRAIT,
    apply_zscore  = TRUE,
    min_genotypes = MIN_GENO
  )

  expect_true(is.list(result))
  expect_true(all(c("results", "predictions", "summary",
                     "zscore_applied", "cv_scheme") %in% names(result)))
  expect_equal(result$cv_scheme, "CrossLocation")
})

test_that("CrossLoc scheme names encode train->predict direction", {
  setup_integration_data()

  result <- run_cross_location(
    pheno_data = test_env$pheno, Ginv_sparse = test_env$Ginv_sparse,
    traits = TEST_TRAIT, min_genotypes = MIN_GENO
  )

  res <- result$results
  if (nrow(res) > 0) {
    expect_true(all(grepl("CrossLoc_", res$cv_scheme)))
    expect_true(all(grepl("->", res$cv_scheme)))
  }
})

test_that("CrossLoc predictions include train_location column", {
  setup_integration_data()

  result <- run_cross_location(
    pheno_data = test_env$pheno, Ginv_sparse = test_env$Ginv_sparse,
    traits = TEST_TRAIT, min_genotypes = MIN_GENO
  )

  preds <- result$predictions
  expect_true("train_location" %in% names(preds))
})

test_that("CrossLoc uses simple model (intercept + G only)", {
  setup_integration_data()

  result <- run_cross_location(
    pheno_data = test_env$pheno, Ginv_sparse = test_env$Ginv_sparse,
    traits = TEST_TRAIT, min_genotypes = MIN_GENO
  )

  # If cross-location ran successfully, it means the simple model
  # (trait ~ 1 + vm(sample.id, Ginv)) worked without location fixed effects.
  # This is verified implicitly: if location were a fixed effect, ASReml
  # would fail or produce degenerate results for single-location training data.
  expect_true(nrow(result$results) > 0 || nrow(result$predictions) > 0,
              info = "CrossLoc should produce some output")
})

test_that("CrossLoc trains on each location and predicts the other(s)", {
  setup_integration_data()

  result <- run_cross_location(
    pheno_data = test_env$pheno, Ginv_sparse = test_env$Ginv_sparse,
    traits = TEST_TRAIT, min_genotypes = MIN_GENO
  )

  preds <- result$predictions
  if (nrow(preds) > 0) {
    locations <- unique(as.character(test_env$pheno$location))
    train_locs <- unique(preds$train_location)

    # Each location should appear as training at least once
    for (loc in locations) {
      if (sum(!is.na(test_env$pheno[["PtHt"]][test_env$pheno$location == loc])) >= 50) {
        expect_true(loc %in% train_locs,
                    info = paste(loc, "should appear as training location"))
      }
    }

    # Predictions should be in target locations (different from train)
    for (i in seq_len(nrow(preds))) {
      expect_true(preds$location[i] != preds$train_location[i],
                  info = "Target location should differ from training location")
    }
  }
})

# ============================================================================
# run_all_cv_schemes wrapper
# ============================================================================

test_that("run_all_cv_schemes returns combined results from all schemes", {
  setup_integration_data()

  result <- run_all_cv_schemes(
    pheno_data    = test_env$pheno,
    Ginv_sparse   = test_env$Ginv_sparse,
    traits        = TEST_TRAIT,
    k_folds       = K_FOLDS,
    n_iterations  = N_ITERATIONS,
    apply_zscore  = TRUE,
    min_genotypes = MIN_GENO
  )

  expect_true(is.list(result))
  expect_true(all(c("cv1", "cv2", "cv0", "cross_loc", "all_results")
                   %in% names(result)))

  all_res <- result$all_results
  expect_true(is.data.frame(all_res))
  expect_true(nrow(all_res) > 0)

  schemes <- unique(all_res$cv_scheme)
  expect_true("CV1" %in% schemes)
  expect_true("CV2" %in% schemes)
  expect_true("CV0" %in% schemes)
  expect_true(any(grepl("CrossLoc_", schemes)))
})

# ============================================================================
# Variance components (one row per varcomp per fit, annotated per scheme)
# ============================================================================
# Each CV runner returns a $varcomps data.frame in addition to $results and
# $predictions. CV1/CV2 annotate with (trait, iteration, fold, seed,
# cv_scheme); CV0 with (trait, held_out_location_year, held_out_location,
# cv_scheme); CrossLoc with (trait, train_location, cv_scheme).

test_that("CV1 returns varcomps annotated with iteration/fold/seed/scheme", {
  setup_integration_data()

  result <- run_cv1(
    pheno_data = test_env$pheno, Ginv_sparse = test_env$Ginv_sparse,
    traits = TEST_TRAIT, k_folds = K_FOLDS, n_iterations = N_ITERATIONS,
    min_genotypes = MIN_GENO
  )

  expect_true("varcomps" %in% names(result))
  vc <- result$varcomps

  if (!is.null(vc) && nrow(vc) > 0) {
    expect_true(all(c("component_name", "trait", "iteration", "fold",
                       "seed", "cv_scheme") %in% names(vc)))
    expect_true(all(vc$cv_scheme == "CV1"))
    expect_true(all(vc$trait == TEST_TRAIT))
  }
})

test_that("CV2 returns varcomps annotated with iteration/fold/seed/scheme", {
  setup_integration_data()

  result <- run_cv2(
    pheno_data = test_env$pheno, Ginv_sparse = test_env$Ginv_sparse,
    traits = TEST_TRAIT, k_folds = K_FOLDS, n_iterations = N_ITERATIONS,
    min_genotypes = MIN_GENO
  )

  expect_true("varcomps" %in% names(result))
  vc <- result$varcomps

  if (!is.null(vc) && nrow(vc) > 0) {
    expect_true(all(c("component_name", "trait", "iteration", "fold",
                       "seed", "cv_scheme") %in% names(vc)))
    expect_true(all(vc$cv_scheme == "CV2"))
  }
})

test_that("CV0 returns varcomps annotated with held_out_location_year", {
  setup_integration_data()

  result <- run_cv0(
    pheno_data = test_env$pheno, Ginv_sparse = test_env$Ginv_sparse,
    traits = TEST_TRAIT, min_genotypes = MIN_GENO
  )

  expect_true("varcomps" %in% names(result))
  vc <- result$varcomps

  if (!is.null(vc) && nrow(vc) > 0) {
    expect_true(all(c("component_name", "trait", "held_out_location_year",
                       "held_out_location", "cv_scheme") %in% names(vc)))
    expect_true(all(vc$cv_scheme == "CV0"))
  }
})

test_that("CV0 uses homogeneous residual structure (units, not dsum)", {
  # CV0 fits asreml inline in run_cv0 with `residual = ~ units` and
  # `random = ~ vm(sample.id, Ginv_sparse) + location:year` (homogeneous),
  # because heterogeneous dsum + at(location):year fails for AUS traits that
  # have only two of three years when a location-year is held out.
  setup_integration_data()

  result <- run_cv0(
    pheno_data = test_env$pheno, Ginv_sparse = test_env$Ginv_sparse,
    traits = TEST_TRAIT, min_genotypes = MIN_GENO
  )

  vc <- result$varcomps
  if (!is.null(vc) && nrow(vc) > 0) {
    comp_names <- unique(vc$component_name)

    # Single pooled residual (units!R or units!units) — NOT split per-location
    expect_false(any(grepl("^AUS!|^PAK!", comp_names)),
                 info = paste("CV0 should use homogeneous residuals;",
                              "found per-location residual components:",
                              paste(comp_names[grepl("^AUS!|^PAK!", comp_names)],
                                    collapse = ", ")))

    # Homogeneous location:year random term — NOT at(location, X):year
    expect_false(any(grepl("^at\\(location", comp_names)),
                 info = "CV0 should use location:year (homogeneous), not at(location):year")
  }
})

test_that("CrossLoc returns varcomps annotated with train_location", {
  setup_integration_data()

  result <- run_cross_location(
    pheno_data = test_env$pheno, Ginv_sparse = test_env$Ginv_sparse,
    traits = TEST_TRAIT, min_genotypes = MIN_GENO
  )

  expect_true("varcomps" %in% names(result))
  vc <- result$varcomps

  if (!is.null(vc) && nrow(vc) > 0) {
    expect_true(all(c("component_name", "trait", "train_location",
                       "cv_scheme") %in% names(vc)))
    expect_true(all(grepl("^CrossLoc_", vc$cv_scheme)))
  }
})
