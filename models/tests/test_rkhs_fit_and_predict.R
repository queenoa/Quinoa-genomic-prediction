library(testthat)

suppressPackageStartupMessages(
  source(file.path("..", "models", "RKHS", "RKHS_utils.R"))
)

# ── Build a minimal but realistic test dataset ────────────────────────────────
#
# Uses the real AUSPAK 1k-marker subset and real phenotype data so we test
# on actual data structure. MCMC is kept very short (100/50/5) for speed.

MARKER_FILE <- file.path("..", "data", "AUSPAK_test_subset_1k.raw")
PHENO_FILE  <- file.path("..", "data", "AUSPAK_phenotypes_means_BLUEs.csv")

# Pre-load data once for all tests (expensive to repeat)
test_dat <- NULL
setup_test_data <- function() {
  if (!is.null(test_dat)) return(test_dat)

  tmpdir <- tempdir()
  kcp <- file.path(tmpdir, "test_fit_kernels.RData")

  dat <- load_and_prepare_data("PtHt_blue", MARKER_FILE,
                                pheno_file = PHENO_FILE,
                                kernel_checkpoint = kcp)

  # Assign to parent environment for caching
  test_dat <<- dat
  dat
}

# Helper: sample n rows evenly across location-years to ensure multiple
# locations are represented (head() only gets the first location)
sample_multi_ly <- function(pheno, n_per_ly = 40) {
  lys <- unique(pheno$location_year)
  sampled <- do.call(rbind, lapply(lys, function(ly) {
    rows <- pheno[pheno$location_year == ly, ]
    rows[seq_len(min(n_per_ly, nrow(rows))), ]
  }))
  sampled
}

# Fast MCMC settings
NITER  <- 100
BURNIN <- 50
THIN   <- 5

# ============================================================================
# CRITICAL TEST: groups argument is NOT compatible with RKHS
# ============================================================================
# BGLR explicitly errors: "model RKHS not implemented for groups"
# fit_rkhs_and_predict has been fixed to not pass groups.

test_that("BGLR RKHS errors when groups argument is passed", {
  dat <- setup_test_data()
  model_data <- sample_multi_ly(dat$pheno, n_per_ly = 40)
  sample_ids <- model_data$sample.id

  K_obs_list <- lapply(dat$K_geno_list, function(K) {
    build_obs_kernel(sample_ids, K)
  })

  Z_fix  <- build_fixed_design(model_data$location, model_data$year)
  groups <- build_groups(model_data$location_year)
  y <- model_data[["PtHt_blue"]]

  ETA <- list(
    list(X = Z_fix,         model = "FIXED"),
    list(K = K_obs_list[[1]], model = "RKHS"),
    list(K = K_obs_list[[2]], model = "RKHS"),
    list(K = K_obs_list[[3]], model = "RKHS")
  )

  # Confirm BGLR rejects groups + RKHS
  expect_error(
    BGLR(y = y, ETA = ETA, groups = groups,
         nIter = NITER, burnIn = BURNIN, thin = THIN,
         verbose = FALSE,
         saveAt = paste0(tempdir(), "/test_groups_err_")),
    "not implemented for groups"
  )

  # Clean up any partial files
  bglr_files <- list.files(tempdir(), pattern = "^test_groups_err_",
                           full.names = TRUE)
  if (length(bglr_files) > 0) file.remove(bglr_files)
})

test_that("fit_rkhs_and_predict works without groups (fixed code)", {
  dat <- setup_test_data()
  model_data <- sample_multi_ly(dat$pheno, n_per_ly = 40)

  # Ensure multiple locations
  expect_gt(length(unique(model_data$location)), 1)

  pred <- fit_rkhs_and_predict(
    model_data  = model_data,
    trait       = "PtHt_blue",
    K_geno_list = dat$K_geno_list,
    nIter = NITER, burnIn = BURNIN, thin = THIN,
    saveAt = paste0(tempdir(), "/test_nogroups_")
  )

  expect_false(is.null(pred),
               info = "fit_rkhs_and_predict should succeed without groups")
  expect_true(is.data.frame(pred))
  expect_gt(nrow(pred), 0)
  expect_true(all(!is.na(pred$predicted)),
              info = "All predictions should be non-NA")
})

test_that("BGLR RKHS without groups produces valid predictions", {
  dat <- setup_test_data()
  model_data <- sample_multi_ly(dat$pheno, n_per_ly = 40)
  sample_ids <- model_data$sample.id

  K_obs_list <- lapply(dat$K_geno_list, function(K) {
    build_obs_kernel(sample_ids, K)
  })

  Z_fix  <- build_fixed_design(model_data$location, model_data$year)
  y <- model_data[["PtHt_blue"]]

  ETA <- list(
    list(X = Z_fix,         model = "FIXED"),
    list(K = K_obs_list[[1]], model = "RKHS"),
    list(K = K_obs_list[[2]], model = "RKHS"),
    list(K = K_obs_list[[3]], model = "RKHS")
  )

  fm <- BGLR(y = y, ETA = ETA,
             nIter = NITER, burnIn = BURNIN, thin = THIN,
             verbose = FALSE,
             saveAt = paste0(tempdir(), "/test_nogroups2_"))

  expect_true(all(!is.na(fm$yHat)))
  expect_equal(length(fm$yHat), nrow(model_data))

  # Variance components should be estimated
  varU <- sapply(fm$ETA[2:4], function(x) x$varU)
  expect_true(all(varU >= 0))

  # Clean up
  bglr_files <- list.files(tempdir(), pattern = "^test_nogroups2_",
                           full.names = TRUE)
  if (length(bglr_files) > 0) file.remove(bglr_files)
})

# ============================================================================
# Predictions structure and validity
# ============================================================================

test_that("fit_rkhs_and_predict returns correct columns", {
  dat <- setup_test_data()
  model_data <- sample_multi_ly(dat$pheno, n_per_ly = 40)

  pred <- fit_rkhs_and_predict(
    model_data  = model_data,
    trait       = "PtHt_blue",
    K_geno_list = dat$K_geno_list,
    nIter = NITER, burnIn = BURNIN, thin = THIN,
    saveAt = paste0(tempdir(), "/test_cols_")
  )

  expect_false(is.null(pred))
  expect_true(all(c("sample.id", "location_year", "location", "predicted")
                   %in% names(pred)))
  expect_equal(nrow(pred), nrow(model_data))
})

test_that("predictions exist for both training and test (NA) rows", {
  dat <- setup_test_data()
  model_data <- sample_multi_ly(dat$pheno, n_per_ly = 40)
  n_rows <- nrow(model_data)

  # Mask some rows
  set.seed(42)
  mask_rows <- sample(1:n_rows, min(40, n_rows %/% 3))
  model_data[["PtHt_blue"]][mask_rows] <- NA

  pred <- fit_rkhs_and_predict(
    model_data  = model_data,
    trait       = "PtHt_blue",
    K_geno_list = dat$K_geno_list,
    nIter = NITER, burnIn = BURNIN, thin = THIN,
    saveAt = paste0(tempdir(), "/test_mask_")
  )

  expect_false(is.null(pred))

  # BGLR predicts ALL rows (including NA/test rows)
  expect_equal(nrow(pred), n_rows)
  expect_true(all(!is.na(pred$predicted)),
              info = "BGLR should predict values for both train and test rows")
})

# ============================================================================
# Kernel variance components
# ============================================================================

test_that("kernel variance components are non-negative for real data", {
  dat <- setup_test_data()
  model_data <- sample_multi_ly(dat$pheno, n_per_ly = 40)
  sample_ids <- model_data$sample.id

  K_obs_list <- lapply(dat$K_geno_list, function(K) {
    build_obs_kernel(sample_ids, K)
  })

  Z_fix <- build_fixed_design(model_data$location, model_data$year)
  y <- model_data[["PtHt_blue"]]

  ETA <- list(
    list(X = Z_fix,         model = "FIXED"),
    list(K = K_obs_list[[1]], model = "RKHS"),
    list(K = K_obs_list[[2]], model = "RKHS"),
    list(K = K_obs_list[[3]], model = "RKHS")
  )

  fm <- BGLR(y = y, ETA = ETA,
             nIter = NITER, burnIn = BURNIN, thin = THIN,
             verbose = FALSE,
             saveAt = paste0(tempdir(), "/test_varu_"))

  varU <- sapply(fm$ETA[2:4], function(x) x$varU)
  expect_true(all(varU >= 0),
              info = "Variance components should be non-negative")

  # Clean up
  bglr_files <- list.files(tempdir(), pattern = "^test_varu_",
                           full.names = TRUE)
  if (length(bglr_files) > 0) file.remove(bglr_files)
})

# ============================================================================
# Error handling
# ============================================================================

test_that("fit_rkhs_and_predict returns NULL when training data too small", {
  dat <- setup_test_data()
  model_data <- sample_multi_ly(dat$pheno, n_per_ly = 5)

  # Set most values to NA so < 50 training observations remain
  n <- nrow(model_data)
  model_data[["PtHt_blue"]][1:(n - 5)] <- NA

  expect_warning(
    pred <- fit_rkhs_and_predict(
      model_data  = model_data,
      trait       = "PtHt_blue",
      K_geno_list = dat$K_geno_list,
      nIter = NITER, burnIn = BURNIN, thin = THIN,
      saveAt = paste0(tempdir(), "/test_small_")
    ),
    "Insufficient training data"
  )

  expect_null(pred)
})

test_that("BGLR temp files are cleaned up after fitting", {
  dat <- setup_test_data()
  model_data <- sample_multi_ly(dat$pheno, n_per_ly = 30)

  # fit_rkhs_and_predict cleans up files matching saveAt prefix in the CWD,
  # so we must use a CWD-relative prefix (matching how the CV scripts work)
  saveAt <- "test_cleanup_rkhs_"

  pred <- fit_rkhs_and_predict(
    model_data  = model_data,
    trait       = "PtHt_blue",
    K_geno_list = dat$K_geno_list,
    nIter = NITER, burnIn = BURNIN, thin = THIN,
    saveAt = saveAt
  )

  expect_false(is.null(pred))

  # Check no BGLR temp files remain in CWD
  remaining <- list.files(pattern = "^test_cleanup_rkhs_", full.names = TRUE)
  expect_equal(length(remaining), 0,
               info = "BGLR temp files should be cleaned up")
})
