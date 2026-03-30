library(testthat)

# Source from GBLUP_utils.R to verify the shared copies are correct
suppressPackageStartupMessages(
  source(file.path("..", "models", "GBLUP", "GBLUP_utils.R"))
)

# ============================================================================
# VALID_TRAITS and LOWER_IS_BETTER_TRAITS constants
# ============================================================================

test_that("VALID_TRAITS matches BayesC/RKHS", {
  expected <- c('DTF_blue', 'DTH_blue', 'PtHt_blue', 'PcleLng_blue',
                'SdLen_blue', 'TGW_blue', 'SdW_z_blue')
  expect_equal(VALID_TRAITS, expected)
})

test_that("LOWER_IS_BETTER_TRAITS matches BayesC/RKHS", {
  expected <- c('DTF_blue', 'DTH_blue', 'PtHt_blue')
  expect_equal(LOWER_IS_BETTER_TRAITS, expected)
})

# ============================================================================
# apply_location_year_scaling
# ============================================================================

test_that("z-scores have mean 0 and sd 1 per location-year", {
  set.seed(42)
  pheno <- data.frame(
    sample.id     = paste0("G", 1:100),
    location_year = rep(c("AUS_2019", "PAK_2020"), each = 50),
    trait1        = c(rnorm(50, mean = 50, sd = 5),
                      rnorm(50, mean = 80, sd = 10)),
    stringsAsFactors = FALSE
  )

  scaled <- apply_location_year_scaling(pheno, "trait1")

  for (ly in c("AUS_2019", "PAK_2020")) {
    vals <- scaled$trait1[scaled$location_year == ly]
    expect_equal(mean(vals), 0, tolerance = 1e-10, label = paste(ly, "mean"))
    expect_equal(sd(vals), 1, tolerance = 1e-10, label = paste(ly, "sd"))
  }
})

test_that("NA values are preserved after scaling", {
  pheno <- data.frame(
    sample.id     = paste0("G", 1:60),
    location_year = rep("AUS_2019", 60),
    trait1        = c(1:55, rep(NA, 5)),
    stringsAsFactors = FALSE
  )

  scaled <- apply_location_year_scaling(pheno, "trait1")
  expect_true(all(is.na(scaled$trait1[56:60])))
})

test_that("scaling warns on too few observations", {
  pheno <- data.frame(
    sample.id     = paste0("G", 1:30),
    location_year = rep("AUS_2019", 30),
    trait1        = c(1:30),
    stringsAsFactors = FALSE
  )

  expect_warning(
    apply_location_year_scaling(pheno, "trait1"),
    "Too few observations"
  )
})

test_that("zero-variance location-year gets sd = 1 (no division by zero)", {
  pheno <- data.frame(
    sample.id     = paste0("G", 1:60),
    location_year = rep("AUS_2019", 60),
    trait1        = rep(5, 60),
    stringsAsFactors = FALSE
  )

  scaled <- apply_location_year_scaling(pheno, "trait1")
  # All values should be 0 (constant - constant) / 1
  expect_true(all(scaled$trait1 == 0))
})

# ============================================================================
# calculate_ndcg
# ============================================================================

test_that("perfect ranking returns NDCG = 1", {
  y_true <- c(10, 8, 6, 4, 2)
  y_pred <- c(10, 8, 6, 4, 2)
  expect_equal(calculate_ndcg(y_true, y_pred, k = 5), 1.0)
})

test_that("lower_is_better flips the ranking", {
  y_true <- c(1, 2, 3, 4, 5)
  y_pred <- c(1, 2, 3, 4, 5)
  result <- calculate_ndcg(y_true, y_pred, k = 5, lower_is_better = TRUE)
  expect_equal(result, 1.0)
})

test_that("reversed ranking gives NDCG < 1", {
  y_true <- c(10, 8, 6, 4, 2)
  y_pred <- c(2, 4, 6, 8, 10)
  result <- calculate_ndcg(y_true, y_pred, k = 5)
  expect_lt(result, 1.0)
  expect_gt(result, 0)
})

test_that("k is clamped to vector length", {
  y_true <- c(5, 3, 1)
  y_pred <- c(5, 3, 1)
  result <- calculate_ndcg(y_true, y_pred, k = 100)
  expect_equal(result, 1.0)
})

# ============================================================================
# evaluate_predictions
# ============================================================================

test_that("evaluate_predictions returns pearson, spearman, ndcg_at_10", {
  set.seed(42)
  y_true <- rnorm(50)
  y_pred <- y_true + rnorm(50, sd = 0.3)

  res <- evaluate_predictions(y_true, y_pred, trait_name = "DTF_blue")

  expect_true(is.list(res))
  expect_named(res, c("pearson", "spearman", "ndcg_at_10"))
  expect_true(res$pearson >= -1 && res$pearson <= 1)
  expect_true(res$spearman >= -1 && res$spearman <= 1)
  expect_true(res$ndcg_at_10 >= 0 && res$ndcg_at_10 <= 1)
})

test_that("evaluate_predictions returns NA for < 2 observations", {
  res <- evaluate_predictions(c(1), c(1))
  expect_true(is.na(res$pearson))
  expect_true(is.na(res$spearman))
  expect_true(is.na(res$ndcg_at_10))
})

test_that("evaluate_predictions uses lower_is_better for NDCG", {
  y_true <- c(1, 2, 3, 4, 5)
  y_pred <- c(1, 2, 3, 4, 5)

  # DTF_blue is in LOWER_IS_BETTER_TRAITS
  res_lower <- evaluate_predictions(y_true, y_pred, trait_name = "DTF_blue")
  # SdW_z_blue is NOT in LOWER_IS_BETTER_TRAITS
  res_higher <- evaluate_predictions(y_true, y_pred, trait_name = "SdW_z_blue")

  # Both should have NDCG = 1 since prediction matches truth perfectly
  expect_equal(res_lower$ndcg_at_10, 1.0)
  expect_equal(res_higher$ndcg_at_10, 1.0)
})

# ============================================================================
# summarise_cv_results
# ============================================================================

test_that("summarise_cv_results returns correct columns", {
  cv_data <- data.frame(
    trait = rep("DTF_blue", 4),
    location_year = rep(c("AUS_2019", "AUS_2020"), each = 2),
    location = rep("AUS", 4),
    pearson = c(0.5, 0.6, 0.7, 0.8),
    spearman = c(0.4, 0.5, 0.6, 0.7),
    ndcg_at_10 = c(0.8, 0.85, 0.9, 0.95),
    n_test_genotypes = rep(50, 4),
    stringsAsFactors = FALSE
  )

  result <- summarise_cv_results(cv_data, "CV1")

  expect_true(is.data.frame(result))
  expected_cols <- c("trait", "location", "pearson_mean", "pearson_sd",
                     "pearson_min", "pearson_max",
                     "spearman_mean", "spearman_sd",
                     "ndcg_at_10_mean", "ndcg_at_10_sd",
                     "n_evaluations", "mean_n_test_genotypes")
  expect_true(all(expected_cols %in% names(result)))
})

test_that("summarise_cv_results warns on empty input", {
  empty_df <- data.frame(
    trait = character(), location_year = character(), location = character(),
    pearson = numeric(), spearman = numeric(), ndcg_at_10 = numeric(),
    n_test_genotypes = integer(), stringsAsFactors = FALSE
  )

  expect_warning(
    summarise_cv_results(empty_df, "CV1"),
    "No successful"
  )
})

test_that("summarise_cv_results groups by trait and location", {
  cv_data <- data.frame(
    trait = rep("DTF_blue", 4),
    location_year = c("AUS_2019", "AUS_2020", "PAK_2019", "PAK_2020"),
    location = c("AUS", "AUS", "PAK", "PAK"),
    pearson = c(0.5, 0.6, 0.3, 0.4),
    spearman = c(0.4, 0.5, 0.2, 0.3),
    ndcg_at_10 = c(0.8, 0.9, 0.7, 0.75),
    n_test_genotypes = rep(50, 4),
    stringsAsFactors = FALSE
  )

  result <- summarise_cv_results(cv_data, "CV1")

  expect_equal(nrow(result), 2)  # AUS and PAK
  aus_row <- result[result$location == "AUS", ]
  pak_row <- result[result$location == "PAK", ]
  expect_equal(aus_row$pearson_mean, 0.55)
  expect_equal(pak_row$pearson_mean, 0.35)
  expect_equal(aus_row$n_evaluations, 2)
})

test_that("summarise_cv_results calculates correct statistics", {
  cv_data <- data.frame(
    trait = rep("DTF_blue", 4),
    location_year = paste0("AUS_", 2019:2022),
    location = rep("AUS", 4),
    pearson = c(0.4, 0.6, 0.8, 1.0),
    spearman = c(0.3, 0.5, 0.7, 0.9),
    ndcg_at_10 = c(0.7, 0.8, 0.9, 1.0),
    n_test_genotypes = rep(50, 4),
    stringsAsFactors = FALSE
  )

  result <- summarise_cv_results(cv_data, "CV1")

  expect_equal(result$pearson_mean, mean(c(0.4, 0.6, 0.8, 1.0)))
  expect_equal(result$pearson_sd, sd(c(0.4, 0.6, 0.8, 1.0)))
  expect_equal(result$pearson_min, 0.4)
  expect_equal(result$pearson_max, 1.0)
  expect_equal(result$n_evaluations, 4)
})
