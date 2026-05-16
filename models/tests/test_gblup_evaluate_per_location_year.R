library(testthat)

suppressPackageStartupMessages(
  source(file.path("..", "GBLUP", "GBLUP_utils.R"))
)

# ============================================================================
# GBLUP evaluate_per_location_year uses join-based approach:
#   - test_ids: character vector of genotype IDs (not integer row indices)
#   - pred_values: must have predicted.value column (ASReml format)
#   - Joins on (sample.id, location_year) to match predictions to observations
#
# This differs from RKHS/BayesC which use index-based lookup.
# ============================================================================

# ── Helper: build synthetic observed and predicted data ─────────────────────

make_test_data <- function(n_per_ly = 20, n_ly = 3, seed = 42) {
  set.seed(seed)
  ly_names <- c("AUS_2019", "PAK_2019", "PAK_2020")[1:n_ly]
  loc_names <- sub("_.*", "", ly_names)

  n <- n_per_ly * n_ly
  geno_ids <- paste0("G", 1:n_per_ly)

  observed_data <- data.frame(
    sample.id     = rep(geno_ids, times = n_ly),
    location_year = rep(ly_names, each = n_per_ly),
    location      = rep(loc_names, each = n_per_ly),
    DTF      = rnorm(n),
    stringsAsFactors = FALSE
  )

  # Simulate ASReml predict output (predicted.value column)
  pred_values <- data.frame(
    sample.id       = rep(geno_ids, times = n_ly),
    location_year   = rep(ly_names, each = n_per_ly),
    predicted.value = observed_data$DTF + rnorm(n, sd = 0.5),
    stringsAsFactors = FALSE
  )

  list(observed = observed_data, predicted = pred_values, geno_ids = geno_ids)
}

# ============================================================================
# Basic functionality
# ============================================================================

test_that("returns list with metrics and predictions", {
  dat <- make_test_data()

  result <- evaluate_per_location_year(
    observed_data = dat$observed,
    pred_values   = dat$predicted,
    test_ids      = dat$geno_ids,
    trait         = "DTF",
    min_genotypes = 5
  )

  expect_true(is.list(result))
  expect_named(result, c("metrics", "predictions"))
  expect_true(is.data.frame(result$metrics))
  expect_true(is.data.frame(result$predictions))
})

test_that("metrics has expected columns", {
  dat <- make_test_data()

  result <- evaluate_per_location_year(
    observed_data = dat$observed,
    pred_values   = dat$predicted,
    test_ids      = dat$geno_ids,
    trait         = "DTF",
    min_genotypes = 5
  )

  expected_cols <- c("trait", "location_year", "location", "pearson",
                     "spearman", "ndcg_at_10", "n_test_genotypes")
  expect_true(all(expected_cols %in% names(result$metrics)))
})

test_that("evaluates all location-years with sufficient data", {
  dat <- make_test_data(n_per_ly = 20, n_ly = 3)

  result <- evaluate_per_location_year(
    observed_data = dat$observed,
    pred_values   = dat$predicted,
    test_ids      = dat$geno_ids,
    trait         = "DTF",
    min_genotypes = 5
  )

  expect_equal(nrow(result$metrics), 3)  # 3 location-years
})

test_that("predictions has expected columns", {
  dat <- make_test_data()

  result <- evaluate_per_location_year(
    observed_data = dat$observed,
    pred_values   = dat$predicted,
    test_ids      = dat$geno_ids,
    trait         = "DTF",
    min_genotypes = 5
  )

  expect_true(all(c("sample.id", "location_year", "location",
                     "observed", "predicted") %in% names(result$predictions)))
})

# ============================================================================
# Join-based matching (key GBLUP difference from RKHS)
# ============================================================================

test_that("only test_ids are included in results", {
  dat <- make_test_data(n_per_ly = 20)
  test_subset <- paste0("G", 1:10)  # Only first 10 genotypes

  result <- evaluate_per_location_year(
    observed_data = dat$observed,
    pred_values   = dat$predicted,
    test_ids      = test_subset,
    trait         = "DTF",
    min_genotypes = 5
  )

  # Predictions should only contain test genotypes
  expect_true(all(result$predictions$sample.id %in% test_subset))
  # Should NOT contain non-test genotypes
  expect_false(any(result$predictions$sample.id %in% paste0("G", 11:20)))
})

test_that("join matches correctly on sample.id AND location_year", {
  set.seed(42)
  # Create data where G1 has different values in different location-years
  observed_data <- data.frame(
    sample.id     = c("G1", "G1", "G2", "G2"),
    location_year = c("AUS_2019", "PAK_2020", "AUS_2019", "PAK_2020"),
    location      = c("AUS", "PAK", "AUS", "PAK"),
    DTF      = c(1.0, 2.0, 3.0, 4.0),
    stringsAsFactors = FALSE
  )

  pred_values <- data.frame(
    sample.id       = c("G1", "G1", "G2", "G2"),
    location_year   = c("AUS_2019", "PAK_2020", "AUS_2019", "PAK_2020"),
    predicted.value = c(1.1, 2.1, 3.1, 4.1),
    stringsAsFactors = FALSE
  )

  result <- evaluate_per_location_year(
    observed_data, pred_values,
    test_ids = c("G1", "G2"), trait = "DTF", min_genotypes = 2
  )

  preds <- result$predictions

  # Verify G1 in AUS_2019 got the right predicted value
  g1_aus <- preds[preds$sample.id == "G1" & preds$location_year == "AUS_2019", ]
  expect_equal(g1_aus$observed, 1.0)
  expect_equal(g1_aus$predicted, 1.1)

  # Verify G1 in PAK_2020 got its own predicted value (not AUS one)
  g1_pak <- preds[preds$sample.id == "G1" & preds$location_year == "PAK_2020", ]
  expect_equal(g1_pak$observed, 2.0)
  expect_equal(g1_pak$predicted, 2.1)
})

# ============================================================================
# min_genotypes filtering
# ============================================================================

test_that("min_genotypes filter skips small location-years", {
  dat <- make_test_data(n_per_ly = 5, n_ly = 1)  # Only 5 genotypes

  result <- evaluate_per_location_year(
    observed_data = dat$observed,
    pred_values   = dat$predicted,
    test_ids      = dat$geno_ids,
    trait         = "DTF",
    min_genotypes = 10  # Higher than available
  )

  expect_equal(nrow(result$metrics), 0)
})

test_that("min_genotypes counts unique genotypes, not rows", {
  # Genotype G1 appears twice in same location-year (shouldn't count double)
  observed_data <- data.frame(
    sample.id     = c("G1", "G1", "G2", "G3"),
    location_year = rep("AUS_2019", 4),
    location      = rep("AUS", 4),
    DTF      = c(1, 1.1, 2, 3),
    stringsAsFactors = FALSE
  )

  pred_values <- data.frame(
    sample.id       = c("G1", "G1", "G2", "G3"),
    location_year   = rep("AUS_2019", 4),
    predicted.value = c(1.1, 1.2, 2.1, 3.1),
    stringsAsFactors = FALSE
  )

  result <- evaluate_per_location_year(
    observed_data, pred_values,
    test_ids = c("G1", "G2", "G3"), trait = "DTF", min_genotypes = 4
  )

  # Only 3 unique genotypes, so min_genotypes = 4 should skip
  expect_equal(nrow(result$metrics), 0)
})

# ============================================================================
# Edge cases
# ============================================================================

test_that("returns empty when no test IDs match predictions", {
  dat <- make_test_data()

  result <- evaluate_per_location_year(
    observed_data = dat$observed,
    pred_values   = dat$predicted,
    test_ids      = c("NONEXISTENT1", "NONEXISTENT2"),
    trait         = "DTF",
    min_genotypes = 5
  )

  expect_equal(nrow(result$metrics), 0)
  expect_equal(nrow(result$predictions), 0)
})

test_that("handles NA trait values in observed data (excludes them)", {
  set.seed(42)
  observed_data <- data.frame(
    sample.id     = paste0("G", 1:20),
    location_year = rep("AUS_2019", 20),
    location      = rep("AUS", 20),
    DTF      = c(rnorm(15), rep(NA, 5)),
    stringsAsFactors = FALSE
  )

  pred_values <- data.frame(
    sample.id       = paste0("G", 1:20),
    location_year   = rep("AUS_2019", 20),
    predicted.value = rnorm(20),
    stringsAsFactors = FALSE
  )

  result <- evaluate_per_location_year(
    observed_data, pred_values,
    test_ids = paste0("G", 1:20), trait = "DTF", min_genotypes = 5
  )

  # Predictions should only include non-NA observations
  expect_equal(nrow(result$predictions), 15)
})

test_that("handles NA predicted values (excludes them from join)", {
  set.seed(42)
  observed_data <- data.frame(
    sample.id     = paste0("G", 1:20),
    location_year = rep("AUS_2019", 20),
    location      = rep("AUS", 20),
    DTF      = rnorm(20),
    stringsAsFactors = FALSE
  )

  pred_values <- data.frame(
    sample.id       = paste0("G", 1:20),
    location_year   = rep("AUS_2019", 20),
    predicted.value = c(rnorm(15), rep(NA, 5)),
    stringsAsFactors = FALSE
  )

  result <- evaluate_per_location_year(
    observed_data, pred_values,
    test_ids = paste0("G", 1:20), trait = "DTF", min_genotypes = 5
  )

  # NA predictions should be filtered out
  expect_equal(nrow(result$predictions), 15)
})

# ============================================================================
# Metrics correctness
# ============================================================================

test_that("metrics match manual calculation", {
  set.seed(42)
  n <- 30
  observed_data <- data.frame(
    sample.id     = paste0("G", 1:n),
    location_year = rep("AUS_2019", n),
    location      = rep("AUS", n),
    DTF      = rnorm(n),
    stringsAsFactors = FALSE
  )

  pred_values <- data.frame(
    sample.id       = paste0("G", 1:n),
    location_year   = rep("AUS_2019", n),
    predicted.value = observed_data$DTF + rnorm(n, sd = 0.3),
    stringsAsFactors = FALSE
  )

  result <- evaluate_per_location_year(
    observed_data, pred_values,
    test_ids = paste0("G", 1:n), trait = "DTF", min_genotypes = 5
  )

  # Manual calculation
  expected_pearson <- cor(observed_data$DTF, pred_values$predicted.value)
  expected_spearman <- cor(observed_data$DTF, pred_values$predicted.value,
                           method = "spearman")

  expect_equal(result$metrics$pearson, expected_pearson, tolerance = 1e-10)
  expect_equal(result$metrics$spearman, expected_spearman, tolerance = 1e-10)
})

test_that("trait column is recorded correctly in metrics", {
  dat <- make_test_data(n_per_ly = 20, n_ly = 1)

  result <- evaluate_per_location_year(
    observed_data = dat$observed,
    pred_values   = dat$predicted,
    test_ids      = dat$geno_ids,
    trait         = "DTF",
    min_genotypes = 5
  )

  expect_true(all(result$metrics$trait == "DTF"))
})

test_that("location is extracted correctly per location-year", {
  dat <- make_test_data(n_per_ly = 20, n_ly = 3)

  result <- evaluate_per_location_year(
    observed_data = dat$observed,
    pred_values   = dat$predicted,
    test_ids      = dat$geno_ids,
    trait         = "DTF",
    min_genotypes = 5
  )

  # AUS_2019 -> AUS, PAK_2019 -> PAK, PAK_2020 -> PAK
  aus_rows <- result$metrics[result$metrics$location_year == "AUS_2019", ]
  pak_rows <- result$metrics[result$metrics$location_year == "PAK_2019", ]
  expect_equal(aus_rows$location, "AUS")
  expect_equal(pak_rows$location, "PAK")
})
