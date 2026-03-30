library(testthat)

suppressPackageStartupMessages(
  source(file.path("..", "models", "BayesC", "BayesC_utils.R"))
)

# Build a synthetic dataset for these tests
make_test_data <- function() {
  n_per_ly <- 20
  observed_data <- data.frame(
    sample.id     = paste0("G", 1:(n_per_ly * 3)),
    location_year = rep(c("AUS_2019", "PAK_2019", "PAK_2020"), each = n_per_ly),
    location      = rep(c("AUS", "PAK", "PAK"), each = n_per_ly),
    DTF_blue      = rnorm(n_per_ly * 3),
    stringsAsFactors = FALSE
  )

  # Simulated predictions (as returned by fit_bayesc_and_predict)
  pred_values <- data.frame(
    sample.id     = observed_data$sample.id,
    location_year = observed_data$location_year,
    location      = observed_data$location,
    predicted     = observed_data$DTF_blue + rnorm(n_per_ly * 3, sd = 0.5),
    stringsAsFactors = FALSE
  )

  list(observed_data = observed_data, pred_values = pred_values)
}

test_that("returns list with metrics and predictions", {
  set.seed(42)
  td <- make_test_data()
  test_indices <- 1:20  # AUS_2019 only

  result <- evaluate_per_location_year(
    td$observed_data, td$pred_values, test_indices,
    trait = "DTF_blue", min_genotypes = 5
  )

  expect_true(is.list(result))
  expect_named(result, c("metrics", "predictions"))
  expect_true(is.data.frame(result$metrics))
  expect_true(is.data.frame(result$predictions))
})

test_that("metrics have correct columns", {
  set.seed(42)
  td <- make_test_data()
  test_indices <- 1:60

  result <- evaluate_per_location_year(
    td$observed_data, td$pred_values, test_indices,
    trait = "DTF_blue", min_genotypes = 5
  )

  expected_cols <- c("trait", "location_year", "location",
                     "pearson", "spearman", "ndcg_at_10", "n_test_genotypes")
  expect_true(all(expected_cols %in% names(result$metrics)))
})

test_that("one row per location-year in metrics", {
  set.seed(42)
  td <- make_test_data()
  test_indices <- 1:60  # all three location-years

  result <- evaluate_per_location_year(
    td$observed_data, td$pred_values, test_indices,
    trait = "DTF_blue", min_genotypes = 5
  )

  expect_equal(nrow(result$metrics), 3)
  expect_equal(sort(result$metrics$location_year),
               c("AUS_2019", "PAK_2019", "PAK_2020"))
})

test_that("min_genotypes filter skips small location-years", {
  set.seed(42)
  td <- make_test_data()
  # Only test 5 genotypes from AUS_2019
  test_indices <- 1:5

  result <- evaluate_per_location_year(
    td$observed_data, td$pred_values, test_indices,
    trait = "DTF_blue", min_genotypes = 10
  )

  # 5 < 10, so AUS_2019 should be skipped -> empty metrics
  # But predictions should still be returned
  expect_equal(nrow(result$metrics), 0)
})

test_that("predictions contain all valid test observations", {
  set.seed(42)
  td <- make_test_data()
  test_indices <- 1:60

  result <- evaluate_per_location_year(
    td$observed_data, td$pred_values, test_indices,
    trait = "DTF_blue", min_genotypes = 5
  )

  # All 60 observations have non-NA observed and predicted
  expect_equal(nrow(result$predictions), 60)
  expect_true(all(c("sample.id", "location_year", "location",
                     "observed", "predicted") %in% names(result$predictions)))
})

test_that("NA observed values are excluded from predictions", {
  set.seed(42)
  td <- make_test_data()
  # Set some observed values to NA
  td$observed_data$DTF_blue[1:5] <- NA
  test_indices <- 1:20

  result <- evaluate_per_location_year(
    td$observed_data, td$pred_values, test_indices,
    trait = "DTF_blue", min_genotypes = 5
  )

  # 20 test indices - 5 NA = 15 valid predictions
  expect_equal(nrow(result$predictions), 15)
  expect_true(all(!is.na(result$predictions$observed)))
})

test_that("returns empty metrics when all observed values are NA", {
  set.seed(42)
  td <- make_test_data()
  td$observed_data$DTF_blue[1:20] <- NA
  test_indices <- 1:20

  result <- evaluate_per_location_year(
    td$observed_data, td$pred_values, test_indices,
    trait = "DTF_blue", min_genotypes = 5
  )

  # No valid pairs -> returns empty metrics data.frame (not a list)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 0)
})

test_that("correlations are between -1 and 1", {
  set.seed(42)
  td <- make_test_data()
  test_indices <- 1:60

  result <- evaluate_per_location_year(
    td$observed_data, td$pred_values, test_indices,
    trait = "DTF_blue", min_genotypes = 5
  )

  expect_true(all(result$metrics$pearson >= -1 & result$metrics$pearson <= 1))
  expect_true(all(result$metrics$spearman >= -1 & result$metrics$spearman <= 1))
  expect_true(all(result$metrics$ndcg_at_10 >= 0 & result$metrics$ndcg_at_10 <= 1))
})
