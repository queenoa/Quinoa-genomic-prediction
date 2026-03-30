library(testthat)

# Lightweight tests for the CV fold assignment logic used by RKHS, BayesC,
# and GBLUP.  Replicates the exact fold creation code from the CV scripts
# and verifies shuffling, balance, and stratification — no model fitting
# required, so these run in under a second.

# ── Synthetic data ──────────────────────────────────────────────────────────

build_test_pheno <- function(n_genotypes = 100, seed = 42) {
  set.seed(seed)
  genotypes <- paste0("G", sprintf("%03d", 1:n_genotypes))
  location_years <- c("AUS_2017", "AUS_2018", "AUS_2019",
                       "PAK_2019", "PAK_2020", "PAK_2021")

  rows <- list()
  for (g in genotypes) {
    n_ly <- sample(3:6, 1)
    lys <- sample(location_years, n_ly)
    for (ly in lys) {
      parts <- strsplit(ly, "_")[[1]]
      rows <- c(rows, list(data.frame(
        sample.id      = g,
        location       = parts[1],
        year           = parts[2],
        location_year  = ly,
        trait_blue     = rnorm(1),
        stringsAsFactors = FALSE
      )))
    }
  }
  do.call(rbind, rows)
}

# ============================================================================
# CV1: Genotype-level fold assignment
#
# Code under test (identical in RKHS_CV1.R:102-103,
# BayesC_CV1_single_iter.R:73-74, GBLUP.R:79-80):
#
#   set.seed(1000 + iter)
#   fold_assignments <- sample(rep(1:k_folds, length.out = n_genotypes))
#   names(fold_assignments) <- genotypes
# ============================================================================

test_that("CV1 folds differ across iterations (shuffling works)", {
  genotypes   <- paste0("G", sprintf("%03d", 1:100))
  k_folds     <- 5
  n_genotypes <- length(genotypes)

  set.seed(1001)
  folds_iter1 <- sample(rep(1:k_folds, length.out = n_genotypes))
  names(folds_iter1) <- genotypes

  set.seed(1002)
  folds_iter2 <- sample(rep(1:k_folds, length.out = n_genotypes))
  names(folds_iter2) <- genotypes

  expect_false(identical(folds_iter1, folds_iter2),
               info = "Fold assignments should differ between iterations")

  fold1_iter1 <- sort(names(folds_iter1)[folds_iter1 == 1])
  fold1_iter2 <- sort(names(folds_iter2)[folds_iter2 == 1])
  expect_false(identical(fold1_iter1, fold1_iter2),
               info = "Fold 1 should contain different genotypes across iterations")
})

test_that("CV1 same seed produces identical fold assignments (reproducible)", {
  genotypes <- paste0("G", sprintf("%03d", 1:100))
  k_folds   <- 5

  set.seed(1001)
  folds_a <- sample(rep(1:k_folds, length.out = length(genotypes)))
  names(folds_a) <- genotypes

  set.seed(1001)
  folds_b <- sample(rep(1:k_folds, length.out = length(genotypes)))
  names(folds_b) <- genotypes

  expect_identical(folds_a, folds_b)
})

test_that("CV1 every genotype is assigned to exactly one fold", {
  genotypes <- paste0("G", sprintf("%03d", 1:100))
  k_folds   <- 5

  set.seed(1001)
  folds <- sample(rep(1:k_folds, length.out = length(genotypes)))
  names(folds) <- genotypes

  expect_equal(length(folds), length(genotypes))
  expect_true(all(folds %in% 1:k_folds))
  expect_equal(length(unique(names(folds))), length(genotypes))
})

test_that("CV1 folds are balanced", {
  genotypes <- paste0("G", sprintf("%03d", 1:100))
  k_folds   <- 5

  set.seed(1001)
  folds <- sample(rep(1:k_folds, length.out = length(genotypes)))

  fold_sizes    <- table(folds)
  expected_size <- length(genotypes) / k_folds

  expect_true(all(fold_sizes >= floor(expected_size)))
  expect_true(all(fold_sizes <= ceiling(expected_size)))
})

test_that("CV1 folds shuffle across many iterations (not just 2)", {
  genotypes   <- paste0("G", sprintf("%03d", 1:100))
  k_folds     <- 5
  n_genotypes <- length(genotypes)

  fold1_sets <- list()
  for (iter in 1:5) {
    set.seed(1000 + iter)
    folds <- sample(rep(1:k_folds, length.out = n_genotypes))
    names(folds) <- genotypes
    fold1_sets[[iter]] <- sort(names(folds)[folds == 1])
  }

  n_unique <- length(unique(fold1_sets))
  expect_equal(n_unique, 5,
               info = "All 5 iterations should produce different fold 1 compositions")
})

# ============================================================================
# CV2: Observation-level fold assignment stratified by location-year
#
# Code under test (identical in RKHS_CV2.R:111-118,
# BayesC_CV2_single_iter.R:85-94, GBLUP.R:235-242):
#
#   set.seed(2000 + iter)
#   fold_assignments <- integer(nrow(pheno))
#   for (ly in unique(location_years_obs)) {
#     ly_obs_idx <- observed_idx[location_years_obs == ly]
#     n_ly <- length(ly_obs_idx)
#     fold_assignments[ly_obs_idx] <- sample(rep(1:k_folds, length.out = n_ly))
#   }
# ============================================================================

test_that("CV2 folds differ across iterations (shuffling works)", {
  pheno              <- build_test_pheno()
  k_folds            <- 5
  observed_idx       <- which(!is.na(pheno$trait_blue))
  location_years_obs <- pheno$location_year[observed_idx]

  assign_cv2_folds <- function(seed) {
    set.seed(seed)
    folds <- integer(nrow(pheno))
    for (ly in unique(location_years_obs)) {
      ly_obs_idx <- observed_idx[location_years_obs == ly]
      folds[ly_obs_idx] <- sample(rep(1:k_folds, length.out = length(ly_obs_idx)))
    }
    folds
  }

  folds1 <- assign_cv2_folds(2001)
  folds2 <- assign_cv2_folds(2002)

  expect_false(identical(folds1[observed_idx], folds2[observed_idx]),
               info = "CV2 fold assignments should differ between iterations")
})

test_that("CV2 same seed produces identical fold assignments (reproducible)", {
  pheno              <- build_test_pheno()
  k_folds            <- 5
  observed_idx       <- which(!is.na(pheno$trait_blue))
  location_years_obs <- pheno$location_year[observed_idx]

  assign_cv2_folds <- function(seed) {
    set.seed(seed)
    folds <- integer(nrow(pheno))
    for (ly in unique(location_years_obs)) {
      ly_obs_idx <- observed_idx[location_years_obs == ly]
      folds[ly_obs_idx] <- sample(rep(1:k_folds, length.out = length(ly_obs_idx)))
    }
    folds
  }

  expect_identical(assign_cv2_folds(2001), assign_cv2_folds(2001))
})

test_that("CV2 all location-years are represented in every fold", {
  pheno              <- build_test_pheno()
  k_folds            <- 5
  observed_idx       <- which(!is.na(pheno$trait_blue))
  location_years_obs <- pheno$location_year[observed_idx]
  all_lys            <- unique(location_years_obs)

  set.seed(2001)
  folds <- integer(nrow(pheno))
  for (ly in all_lys) {
    ly_obs_idx <- observed_idx[location_years_obs == ly]
    folds[ly_obs_idx] <- sample(rep(1:k_folds, length.out = length(ly_obs_idx)))
  }

  for (fold in 1:k_folds) {
    fold_idx    <- which(folds == fold)
    lys_in_fold <- unique(pheno$location_year[fold_idx])
    expect_true(setequal(lys_in_fold, all_lys),
                info = paste("Fold", fold, "should contain all location-years.",
                             "Missing:", paste(setdiff(all_lys, lys_in_fold),
                                               collapse = ", ")))
  }
})

test_that("CV2 folds are balanced within each location-year", {
  pheno              <- build_test_pheno()
  k_folds            <- 5
  observed_idx       <- which(!is.na(pheno$trait_blue))
  location_years_obs <- pheno$location_year[observed_idx]

  set.seed(2001)
  folds <- integer(nrow(pheno))
  for (ly in unique(location_years_obs)) {
    ly_obs_idx <- observed_idx[location_years_obs == ly]
    folds[ly_obs_idx] <- sample(rep(1:k_folds, length.out = length(ly_obs_idx)))
  }

  for (ly in unique(location_years_obs)) {
    ly_idx     <- which(pheno$location_year == ly & folds > 0)
    fold_sizes <- table(folds[ly_idx])
    n_ly       <- length(ly_idx)
    expected   <- n_ly / k_folds

    expect_true(all(fold_sizes >= floor(expected)),
                info = paste(ly, ": fold too small"))
    expect_true(all(fold_sizes <= ceiling(expected)),
                info = paste(ly, ": fold too large"))
  }
})

test_that("CV2 unobserved rows remain unassigned (fold = 0)", {
  pheno <- build_test_pheno()
  # Inject some NAs
  pheno$trait_blue[sample(nrow(pheno), 20)] <- NA

  k_folds            <- 5
  observed_idx       <- which(!is.na(pheno$trait_blue))
  unobserved_idx     <- which(is.na(pheno$trait_blue))
  location_years_obs <- pheno$location_year[observed_idx]

  set.seed(2001)
  folds <- integer(nrow(pheno))
  for (ly in unique(location_years_obs)) {
    ly_obs_idx <- observed_idx[location_years_obs == ly]
    folds[ly_obs_idx] <- sample(rep(1:k_folds, length.out = length(ly_obs_idx)))
  }

  expect_true(all(folds[unobserved_idx] == 0),
              info = "Unobserved rows should not be assigned to any fold")
  expect_true(all(folds[observed_idx] > 0),
              info = "All observed rows should be assigned to a fold")
})

test_that("CV2 folds shuffle across many iterations (not just 2)", {
  pheno              <- build_test_pheno()
  k_folds            <- 5
  observed_idx       <- which(!is.na(pheno$trait_blue))
  location_years_obs <- pheno$location_year[observed_idx]

  fingerprints <- character(5)
  for (iter in 1:5) {
    set.seed(2000 + iter)
    folds <- integer(nrow(pheno))
    for (ly in unique(location_years_obs)) {
      ly_obs_idx <- observed_idx[location_years_obs == ly]
      folds[ly_obs_idx] <- sample(rep(1:k_folds, length.out = length(ly_obs_idx)))
    }
    fingerprints[iter] <- paste(folds[observed_idx], collapse = ",")
  }

  expect_equal(length(unique(fingerprints)), 5,
               info = "All 5 iterations should produce distinct fold assignments")
})
