library(testthat)

suppressPackageStartupMessages(
  source(file.path("..", "models", "RKHS", "RKHS_utils.R"))
)

# ── Synthetic marker matrix for kernel tests ─────────────────────────────────

make_test_markers <- function(n_geno = 20, n_snps = 50, seed = 42) {
  set.seed(seed)
  X <- matrix(sample(0:2, n_geno * n_snps, replace = TRUE),
              nrow = n_geno, ncol = n_snps)
  rownames(X) <- paste0("G", 1:n_geno)
  colnames(X) <- paste0("SNP", 1:n_snps)
  X
}

# Helper: compute kernels from raw marker matrix (replicating load_and_prepare_data logic)
compute_test_kernels <- function(X) {
  X <- scale(X, center = TRUE, scale = TRUE)

  # Remove monomorphic markers (NaN after scaling)
  nan_cols <- which(colSums(is.nan(X)) > 0)
  if (length(nan_cols) > 0) X <- X[, -nan_cols]

  D <- as.matrix(dist(X, method = "euclidean"))^2
  D <- D / mean(D)

  med_D <- median(D[lower.tri(D)])
  h_base <- 1 / med_D
  h_values <- h_base * c(1/5, 1, 5)

  K_list <- lapply(h_values, function(h) {
    K <- exp(-h * D)
    rownames(K) <- colnames(K) <- rownames(X)
    K
  })

  list(K_list = K_list, h_values = h_values, D = D, med_D = med_D)
}

# ============================================================================
# Distance matrix properties
# ============================================================================

test_that("squared Euclidean distance matrix is symmetric", {
  X <- make_test_markers()
  res <- compute_test_kernels(X)
  expect_true(isSymmetric(res$D))
})

test_that("distance matrix diagonal is zero", {
  X <- make_test_markers()
  res <- compute_test_kernels(X)
  expect_equal(unname(diag(res$D)), rep(0, nrow(X)))
})

test_that("all off-diagonal distances are positive", {
  X <- make_test_markers()
  res <- compute_test_kernels(X)
  off_diag <- res$D[lower.tri(res$D)]
  expect_true(all(off_diag > 0))
})

test_that("normalised distance matrix has mean 1", {
  X <- make_test_markers()
  res <- compute_test_kernels(X)
  expect_equal(mean(res$D), 1, tolerance = 1e-10)
})

# ============================================================================
# Gaussian kernel properties
# ============================================================================

test_that("all three kernels are symmetric", {
  X <- make_test_markers()
  res <- compute_test_kernels(X)
  for (i in 1:3) {
    expect_true(isSymmetric(res$K_list[[i]]),
                info = paste("Kernel", i, "should be symmetric"))
  }
})

test_that("kernel diagonal is exactly 1 (exp(-h * 0) = 1)", {
  X <- make_test_markers()
  res <- compute_test_kernels(X)
  for (i in 1:3) {
    expect_equal(unname(diag(res$K_list[[i]])), rep(1, nrow(X)),
                 info = paste("Kernel", i, "diagonal should be 1"))
  }
})

test_that("all kernel values are in (0, 1]", {
  X <- make_test_markers()
  res <- compute_test_kernels(X)
  for (i in 1:3) {
    K <- res$K_list[[i]]
    expect_true(all(K > 0), info = paste("Kernel", i, ": all values > 0"))
    expect_true(all(K <= 1), info = paste("Kernel", i, ": all values <= 1"))
  }
})

test_that("kernels are positive semi-definite", {
  X <- make_test_markers()
  res <- compute_test_kernels(X)
  for (i in 1:3) {
    eigenvalues <- eigen(res$K_list[[i]], symmetric = TRUE, only.values = TRUE)$values
    # Allow small numerical tolerance for near-zero eigenvalues
    expect_true(all(eigenvalues > -1e-10),
                info = paste("Kernel", i, "should be PSD"))
  }
})

test_that("genotype IDs are preserved as row/colnames", {
  X <- make_test_markers(n_geno = 15)
  res <- compute_test_kernels(X)
  expected_ids <- paste0("G", 1:15)
  for (i in 1:3) {
    expect_equal(rownames(res$K_list[[i]]), expected_ids)
    expect_equal(colnames(res$K_list[[i]]), expected_ids)
  }
})

# ============================================================================
# Bandwidth ordering
# ============================================================================

test_that("three bandwidths span 25-fold range with correct ratios", {
  X <- make_test_markers()
  res <- compute_test_kernels(X)

  h <- res$h_values
  expect_equal(length(h), 3)
  # h1 < h2 < h3

expect_lt(h[1], h[2])
  expect_lt(h[2], h[3])
  # h3/h1 = 25
  expect_equal(h[3] / h[1], 25, tolerance = 1e-10)
  # h2/h1 = 5
  expect_equal(h[2] / h[1], 5, tolerance = 1e-10)
})

test_that("wider bandwidth (larger h) produces smaller off-diagonal values", {
  # Larger h -> more decay -> values closer to 0
  X <- make_test_markers()
  res <- compute_test_kernels(X)

  mean_offdiag <- sapply(res$K_list, function(K) mean(K[lower.tri(K)]))
  # K_wide (h1, small h) should have largest off-diagonal
  # K_narrow (h3, large h) should have smallest off-diagonal
  expect_gt(mean_offdiag[1], mean_offdiag[2])
  expect_gt(mean_offdiag[2], mean_offdiag[3])
})

# ============================================================================
# Edge cases
# ============================================================================

test_that("identical genotypes have kernel value 1", {
  # Create matrix where rows 1 and 2 are identical
  X <- make_test_markers(n_geno = 10)
  X[2, ] <- X[1, ]
  res <- compute_test_kernels(X)

  for (i in 1:3) {
    expect_equal(res$K_list[[i]][1, 2], 1.0,
                 info = paste("Kernel", i, ": identical genotypes should have K=1"))
  }
})

test_that("monomorphic markers are removed before kernel computation", {
  X <- make_test_markers(n_geno = 20, n_snps = 50)
  # Make 5 columns monomorphic
  X[, 1:5] <- 1

  # scale() will produce NaN for zero-variance columns
  X_scaled <- scale(X, center = TRUE, scale = TRUE)
  nan_cols <- which(colSums(is.nan(X_scaled)) > 0)
  expect_equal(length(nan_cols), 5)

  # compute_test_kernels handles this
  res <- compute_test_kernels(X)
  expect_equal(nrow(res$K_list[[1]]), 20)
  expect_true(all(!is.nan(res$K_list[[1]])))
})
