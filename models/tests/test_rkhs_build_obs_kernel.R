library(testthat)

suppressPackageStartupMessages(
  source(file.path("..", "RKHS", "RKHS_utils.R"))
)

# ── Helper: make a small genotype-level kernel ────────────────────────────────

make_geno_kernel <- function(n_geno = 5) {
  ids <- paste0("G", 1:n_geno)
  # Simple kernel: identity + small off-diagonal
  K <- diag(n_geno) * 0.5 + 0.5
  rownames(K) <- colnames(K) <- ids
  K
}

# ============================================================================
# Basic functionality
# ============================================================================

test_that("build_obs_kernel returns correct dimensions for unique IDs", {
  K_geno <- make_geno_kernel(5)
  sample_ids <- c("G1", "G2", "G3")
  K_obs <- build_obs_kernel(sample_ids, K_geno)

  expect_equal(nrow(K_obs), 3)
  expect_equal(ncol(K_obs), 3)
  expect_true(is.matrix(K_obs))
})

test_that("build_obs_kernel preserves kernel values", {
  K_geno <- make_geno_kernel(5)
  # Set a distinctive value
  K_geno["G2", "G3"] <- 0.42
  K_geno["G3", "G2"] <- 0.42

  sample_ids <- c("G1", "G2", "G3")
  K_obs <- build_obs_kernel(sample_ids, K_geno)

  expect_equal(K_obs[2, 3], 0.42)
  expect_equal(K_obs[3, 2], 0.42)
})

test_that("build_obs_kernel clears row/colnames", {
  K_geno <- make_geno_kernel(5)
  sample_ids <- c("G1", "G2", "G3")
  K_obs <- build_obs_kernel(sample_ids, K_geno)

  expect_null(rownames(K_obs))
  expect_null(colnames(K_obs))
})

# ============================================================================
# Multi-environment expansion (duplicated genotypes)
# ============================================================================

test_that("duplicated genotypes expand kernel correctly", {
  K_geno <- make_geno_kernel(3)
  # Set distinctive values
  K_geno["G1", "G2"] <- 0.7
  K_geno["G2", "G1"] <- 0.7
  K_geno["G1", "G3"] <- 0.3
  K_geno["G3", "G1"] <- 0.3

  # G1 appears in two environments, G2 once, G3 once
  sample_ids <- c("G1", "G2", "G3", "G1")
  K_obs <- build_obs_kernel(sample_ids, K_geno)

  expect_equal(nrow(K_obs), 4)
  expect_equal(ncol(K_obs), 4)

  # K_obs[1,4] and K_obs[4,1] should be K_geno["G1","G1"] = 1.0
  expect_equal(K_obs[1, 4], K_geno["G1", "G1"])
  expect_equal(K_obs[4, 1], K_geno["G1", "G1"])

  # K_obs[1,2] and K_obs[4,2] should both be K_geno["G1","G2"] = 0.7
  expect_equal(K_obs[1, 2], 0.7)
  expect_equal(K_obs[4, 2], 0.7)

  # K_obs[1,3] and K_obs[4,3] should both be K_geno["G1","G3"] = 0.3
  expect_equal(K_obs[1, 3], 0.3)
  expect_equal(K_obs[4, 3], 0.3)
})

test_that("observation-level kernel is symmetric even with duplicates", {
  K_geno <- make_geno_kernel(4)
  # Multi-env: each genotype appears 3 times
  sample_ids <- rep(paste0("G", 1:4), each = 3)
  K_obs <- build_obs_kernel(sample_ids, K_geno)

  expect_equal(nrow(K_obs), 12)
  expect_true(isSymmetric(K_obs))
})

test_that("observation-level kernel diagonal matches genotype self-similarity", {
  n_geno <- 5
  K_geno <- make_geno_kernel(n_geno)

  # Multi-env data: genotypes repeated across 3 location-years
  sample_ids <- rep(paste0("G", 1:n_geno), times = 3)
  K_obs <- build_obs_kernel(sample_ids, K_geno)

  # All diagonal values should equal K_geno[gi,gi]
  for (i in seq_along(sample_ids)) {
    gid <- sample_ids[i]
    expect_equal(K_obs[i, i], K_geno[gid, gid],
                 info = paste("Row", i, "genotype", gid))
  }
})

# ============================================================================
# Full-size test with real-ish data
# ============================================================================

test_that("build_obs_kernel scales to realistic multi-environment data", {
  # 100 genotypes, 6 location-years -> ~600 observations
  n_geno <- 100
  K_geno <- diag(n_geno) * 0.3 + 0.7
  ids <- paste0("G", 1:n_geno)
  rownames(K_geno) <- colnames(K_geno) <- ids

  sample_ids <- rep(ids, times = 6)  # 600 observations
  K_obs <- build_obs_kernel(sample_ids, K_geno)

  expect_equal(nrow(K_obs), 600)
  expect_equal(ncol(K_obs), 600)
  expect_true(is.matrix(K_obs))
  expect_null(rownames(K_obs))
})

# ============================================================================
# Edge cases
# ============================================================================

test_that("single observation returns 1x1 kernel", {
  K_geno <- make_geno_kernel(5)
  K_obs <- build_obs_kernel("G3", K_geno)

  expect_equal(nrow(K_obs), 1)
  expect_equal(ncol(K_obs), 1)
  expect_equal(K_obs[1, 1], K_geno["G3", "G3"])
})

test_that("all same genotype returns matrix of identical values", {
  K_geno <- make_geno_kernel(3)
  sample_ids <- rep("G2", 5)
  K_obs <- build_obs_kernel(sample_ids, K_geno)

  # All values should be K_geno["G2","G2"]
  expected_val <- K_geno["G2", "G2"]
  expect_true(all(K_obs == expected_val))
})
