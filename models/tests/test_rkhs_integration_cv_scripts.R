library(testthat)

# Integration tests for RKHS CV scripts.
# Runs each script end-to-end with 1k-marker test data and minimal MCMC.
# These tests take ~2-5 minutes total due to BGLR fitting.

MODELS_DIR  <- normalizePath(file.path("..", "models", "RKHS"))
MARKER_FILE <- normalizePath(file.path("..", "data", "AUSPAK_test_subset_1k.raw"))
PHENO_FILE  <- normalizePath(file.path("..", "data", "AUSPAK_phenotypes_means_BLUEs.csv"))

# Use a temp directory for all output files
OUTDIR <- file.path(tempdir(), "rkhs_integration_tests")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Environment variables for fast MCMC and single iteration
MCMC_ENV <- c(
  RKHS_NITER    = "100",
  RKHS_BURNIN   = "50",
  RKHS_THIN     = "5",
  RKHS_CV_ITERS = "2"
)

# Test trait — use one with data in most location-years
TEST_TRAIT <- "PtHt_blue"

# Helper: run an Rscript in OUTDIR with fast MCMC settings
run_cv_script <- function(script_name, args) {
  script_path <- file.path(MODELS_DIR, script_name)
  # Copy the utils file to OUTDIR so source("RKHS_utils.R") works
  file.copy(file.path(MODELS_DIR, "RKHS_utils.R"),
            file.path(OUTDIR, "RKHS_utils.R"), overwrite = TRUE)

  cmd <- paste(
    paste(paste0(names(MCMC_ENV), "=", MCMC_ENV), collapse = " "),
    "Rscript", shQuote(script_path),
    paste(shQuote(args), collapse = " ")
  )
  result <- system(cmd, intern = TRUE, ignore.stderr = FALSE)
  attr(result, "status")
}

# ============================================================================
# CV1
# ============================================================================

test_that("CV1 script runs and produces expected output files", {
  old_wd <- setwd(OUTDIR)
  on.exit(setwd(old_wd))

  # Copy phenotype file so default path works
  file.copy(PHENO_FILE, file.path(OUTDIR, "AUSPAK_phenotypes_means_BLUEs.csv"),
            overwrite = TRUE)

  status <- run_cv_script("RKHS_CV1.R", c(TEST_TRAIT, MARKER_FILE))

  results_file <- sprintf("cv_results_CV1_%s_RKHS.csv", TEST_TRAIT)
  preds_file   <- sprintf("predictions_CV1_%s_RKHS.csv", TEST_TRAIT)

  expect_true(file.exists(results_file), info = "CV1 results CSV should exist")
  expect_true(file.exists(preds_file), info = "CV1 predictions CSV should exist")
})

test_that("CV1 results CSV has correct structure", {
  old_wd <- setwd(OUTDIR)
  on.exit(setwd(old_wd))

  results_file <- sprintf("cv_results_CV1_%s_RKHS.csv", TEST_TRAIT)
  res <- read.csv(results_file, stringsAsFactors = FALSE)

  expected_cols <- c("iteration", "fold", "trait", "location_year", "location",
                     "pearson", "spearman", "ndcg_at_10", "seed",
                     "n_test_genotypes", "cv_scheme")
  expect_true(all(expected_cols %in% names(res)),
              info = paste("Missing columns:", paste(setdiff(expected_cols, names(res)),
                                                      collapse = ", ")))

  expect_true(nrow(res) > 0, info = "Should have at least one evaluation row")
  expect_true(all(res$cv_scheme == "CV1"))
  expect_true(all(res$trait == TEST_TRAIT))
  expect_true(all(res$pearson >= -1 & res$pearson <= 1))
  expect_true(all(res$ndcg_at_10 >= 0 & res$ndcg_at_10 <= 1))
})

test_that("CV1 predictions CSV has correct structure", {
  old_wd <- setwd(OUTDIR)
  on.exit(setwd(old_wd))

  preds_file <- sprintf("predictions_CV1_%s_RKHS.csv", TEST_TRAIT)
  preds <- read.csv(preds_file, stringsAsFactors = FALSE)

  expected_cols <- c("sample.id", "location_year", "location",
                     "observed", "predicted", "trait", "iteration",
                     "fold", "seed", "cv_scheme")
  expect_true(all(expected_cols %in% names(preds)))

  expect_true(nrow(preds) > 0)
  expect_true(all(!is.na(preds$observed)))
  expect_true(all(!is.na(preds$predicted)))
  expect_true(all(preds$cv_scheme == "CV1"))
})

test_that("CV1 seed scheme matches BayesC (1000 + iter)", {
  old_wd <- setwd(OUTDIR)
  on.exit(setwd(old_wd))

  results_file <- sprintf("cv_results_CV1_%s_RKHS.csv", TEST_TRAIT)
  res <- read.csv(results_file, stringsAsFactors = FALSE)

  for (iter in unique(res$iteration)) {
    iter_seeds <- unique(res$seed[res$iteration == iter])
    expect_equal(iter_seeds, 1000 + iter,
                 info = paste("Iteration", iter, "seed should be", 1000 + iter))
  }
})

test_that("CV1 fold assignments differ between iterations (shuffling)", {
  old_wd <- setwd(OUTDIR)
  on.exit(setwd(old_wd))

  preds_file <- sprintf("predictions_CV1_%s_RKHS.csv", TEST_TRAIT)
  preds <- read.csv(preds_file, stringsAsFactors = FALSE)

  # Need at least 2 iterations to compare
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
# CV2
# ============================================================================

test_that("CV2 script runs and produces expected output files", {
  old_wd <- setwd(OUTDIR)
  on.exit(setwd(old_wd))

  status <- run_cv_script("RKHS_CV2.R", c(TEST_TRAIT, MARKER_FILE))

  results_file <- sprintf("cv_results_CV2_%s_RKHS.csv", TEST_TRAIT)
  preds_file   <- sprintf("predictions_CV2_%s_RKHS.csv", TEST_TRAIT)

  expect_true(file.exists(results_file), info = "CV2 results CSV should exist")
  expect_true(file.exists(preds_file), info = "CV2 predictions CSV should exist")
})

test_that("CV2 results CSV has correct structure", {
  old_wd <- setwd(OUTDIR)
  on.exit(setwd(old_wd))

  results_file <- sprintf("cv_results_CV2_%s_RKHS.csv", TEST_TRAIT)
  res <- read.csv(results_file, stringsAsFactors = FALSE)

  expect_true(nrow(res) > 0)
  expect_true(all(res$cv_scheme == "CV2"))
  expect_true(all(res$trait == TEST_TRAIT))
  expect_true(all(res$pearson >= -1 & res$pearson <= 1))
})

test_that("CV2 seed scheme matches BayesC (2000 + iter)", {
  old_wd <- setwd(OUTDIR)
  on.exit(setwd(old_wd))

  results_file <- sprintf("cv_results_CV2_%s_RKHS.csv", TEST_TRAIT)
  res <- read.csv(results_file, stringsAsFactors = FALSE)

  for (iter in unique(res$iteration)) {
    iter_seeds <- unique(res$seed[res$iteration == iter])
    expect_equal(iter_seeds, 2000 + iter,
                 info = paste("Iteration", iter, "seed should be", 2000 + iter))
  }
})

test_that("CV2 predictions CSV has correct structure", {
  old_wd <- setwd(OUTDIR)
  on.exit(setwd(old_wd))

  preds_file <- sprintf("predictions_CV2_%s_RKHS.csv", TEST_TRAIT)
  preds <- read.csv(preds_file, stringsAsFactors = FALSE)

  expect_true(nrow(preds) > 0)
  expect_true(all(!is.na(preds$observed)))
  expect_true(all(!is.na(preds$predicted)))
  expect_true(all(preds$cv_scheme == "CV2"))
})

test_that("CV2 all location-years have predictions in every fold", {
  old_wd <- setwd(OUTDIR)
  on.exit(setwd(old_wd))

  preds_file <- sprintf("predictions_CV2_%s_RKHS.csv", TEST_TRAIT)
  preds <- read.csv(preds_file, stringsAsFactors = FALSE)

  # Check within each iteration
  for (iter in unique(preds$iteration)) {
    iter_preds <- preds[preds$iteration == iter, ]
    all_lys    <- unique(iter_preds$location_year)
    k_folds    <- max(iter_preds$fold)

    for (fold in 1:k_folds) {
      fold_lys <- unique(iter_preds$location_year[iter_preds$fold == fold])
      expect_true(setequal(fold_lys, all_lys),
                  info = paste("Iter", iter, "fold", fold,
                               "should have predictions for all location-years.",
                               "Missing:",
                               paste(setdiff(all_lys, fold_lys), collapse = ", ")))
    }
  }
})

test_that("CV2 fold assignments differ between iterations (shuffling)", {
  old_wd <- setwd(OUTDIR)
  on.exit(setwd(old_wd))

  preds_file <- sprintf("predictions_CV2_%s_RKHS.csv", TEST_TRAIT)
  preds <- read.csv(preds_file, stringsAsFactors = FALSE)

  iters <- sort(unique(preds$iteration))
  expect_gte(length(iters), 2,
             label = "Need >= 2 iterations to test fold shuffling")

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
})

# ============================================================================
# CV0 + CrossLoc
# ============================================================================

test_that("CV0+CrossLoc script runs and produces expected output files", {
  old_wd <- setwd(OUTDIR)
  on.exit(setwd(old_wd))

  status <- run_cv_script("RKHS_CV0_CrossLoc.R", c(TEST_TRAIT, MARKER_FILE))

  expect_true(file.exists(sprintf("cv_results_CV0_%s_RKHS.csv", TEST_TRAIT)),
              info = "CV0 results CSV should exist")
  expect_true(file.exists(sprintf("predictions_CV0_%s_RKHS.csv", TEST_TRAIT)),
              info = "CV0 predictions CSV should exist")
  expect_true(file.exists(sprintf("cv_results_CrossLoc_%s_RKHS.csv", TEST_TRAIT)),
              info = "CrossLoc results CSV should exist")
  expect_true(file.exists(sprintf("predictions_CrossLoc_%s_RKHS.csv", TEST_TRAIT)),
              info = "CrossLoc predictions CSV should exist")
})

test_that("CV0 results CSV has correct structure", {
  old_wd <- setwd(OUTDIR)
  on.exit(setwd(old_wd))

  res <- read.csv(sprintf("cv_results_CV0_%s_RKHS.csv", TEST_TRAIT),
                   stringsAsFactors = FALSE)

  expect_true(nrow(res) > 0)
  expect_true(all(res$cv_scheme == "CV0"))
  expect_true(all(res$trait == TEST_TRAIT))
  expect_true(all(res$pearson >= -1 & res$pearson <= 1))

  # CV0 is leave-one-location-year-out: should have results for multiple LYs
  expect_gt(length(unique(res$location_year)), 1)
})

test_that("CV0 predictions CSV has correct structure", {
  old_wd <- setwd(OUTDIR)
  on.exit(setwd(old_wd))

  preds <- read.csv(sprintf("predictions_CV0_%s_RKHS.csv", TEST_TRAIT),
                     stringsAsFactors = FALSE)

  expect_true(nrow(preds) > 0)
  expect_true(all(!is.na(preds$observed)))
  expect_true(all(!is.na(preds$predicted)))
})

test_that("CrossLoc results CSV has correct structure", {
  old_wd <- setwd(OUTDIR)
  on.exit(setwd(old_wd))

  res <- read.csv(sprintf("cv_results_CrossLoc_%s_RKHS.csv", TEST_TRAIT),
                   stringsAsFactors = FALSE)

  expect_true(nrow(res) > 0)
  expect_true(all(res$trait == TEST_TRAIT))
  expect_true(all(res$pearson >= -1 & res$pearson <= 1))

  # CrossLoc cv_scheme encodes train->predict direction
  expect_true(all(grepl("CrossLoc_", res$cv_scheme)))
})

test_that("CrossLoc predictions CSV has correct structure", {
  old_wd <- setwd(OUTDIR)
  on.exit(setwd(old_wd))

  preds <- read.csv(sprintf("predictions_CrossLoc_%s_RKHS.csv", TEST_TRAIT),
                     stringsAsFactors = FALSE)

  expect_true(nrow(preds) > 0)
  expect_true(all(!is.na(preds$observed)))
  expect_true(all(!is.na(preds$predicted)))
  expect_true("train_location" %in% names(preds))
})

# ============================================================================
# Aggregation (consumes CSVs from above)
# ============================================================================

test_that("aggregate script produces combined output", {
  old_wd <- setwd(OUTDIR)
  on.exit(setwd(old_wd))

  agg_script <- file.path(MODELS_DIR, "aggregate_RKHS_results.R")
  cmd <- paste("Rscript", shQuote(agg_script), TEST_TRAIT)
  system(cmd, intern = TRUE)

  combined_file <- sprintf("cv_results_%s_RKHS_all_schemes.csv", TEST_TRAIT)
  summary_file  <- sprintf("cv_summary_%s_RKHS_all_schemes.csv", TEST_TRAIT)
  preds_file    <- sprintf("predictions_%s_RKHS_all_schemes.csv", TEST_TRAIT)

  expect_true(file.exists(combined_file))
  expect_true(file.exists(summary_file))
  expect_true(file.exists(preds_file))
})

test_that("aggregated results contain all CV schemes", {
  old_wd <- setwd(OUTDIR)
  on.exit(setwd(old_wd))

  combined_file <- sprintf("cv_results_%s_RKHS_all_schemes.csv", TEST_TRAIT)
  res <- read.csv(combined_file, stringsAsFactors = FALSE)

  schemes <- unique(res$cv_scheme)
  expect_true("CV1" %in% schemes)
  expect_true("CV2" %in% schemes)
  expect_true("CV0" %in% schemes)
  expect_true(any(grepl("CrossLoc_", schemes)))
})

test_that("summary has expected columns and grouping", {
  old_wd <- setwd(OUTDIR)
  on.exit(setwd(old_wd))

  summary_file <- sprintf("cv_summary_%s_RKHS_all_schemes.csv", TEST_TRAIT)
  summ <- read.csv(summary_file, stringsAsFactors = FALSE)

  expect_true(all(c("trait", "cv_scheme", "location",
                     "pearson_mean", "pearson_sd",
                     "n_evaluations") %in% names(summ)))
  expect_true(nrow(summ) > 0)
  expect_true(all(summ$trait == TEST_TRAIT))
})

# ============================================================================
# RKHS-specific checks: non-empty predictions confirm BGLR actually ran
# ============================================================================

test_that("CV1 predictions have reasonable variance (not constant)", {
  old_wd <- setwd(OUTDIR)
  on.exit(setwd(old_wd))

  preds_file <- sprintf("predictions_CV1_%s_RKHS.csv", TEST_TRAIT)
  preds <- read.csv(preds_file, stringsAsFactors = FALSE)

  # If BGLR silently fails (e.g., groups issue), predictions might be constant
  pred_sd <- sd(preds$predicted)
  expect_gt(pred_sd, 0.01,
            label = "Prediction variance should be non-trivial (not constant)")
})

test_that("CV2 predictions have reasonable variance (not constant)", {
  old_wd <- setwd(OUTDIR)
  on.exit(setwd(old_wd))

  preds_file <- sprintf("predictions_CV2_%s_RKHS.csv", TEST_TRAIT)
  preds <- read.csv(preds_file, stringsAsFactors = FALSE)

  pred_sd <- sd(preds$predicted)
  expect_gt(pred_sd, 0.01,
            label = "Prediction variance should be non-trivial")
})

# Clean up temp directory
test_that("cleanup", {
  unlink(OUTDIR, recursive = TRUE)
  expect_true(TRUE)
})
