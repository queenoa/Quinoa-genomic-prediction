library(testthat)

suppressPackageStartupMessages(
  source(file.path("..", "models", "BayesC", "BayesC_utils.R"))
)

test_that("returns correct Pearson and Spearman for perfect correlation", {
  y_true <- c(1, 2, 3, 4, 5)
  y_pred <- c(1, 2, 3, 4, 5)
  res <- evaluate_predictions(y_true, y_pred)
  expect_equal(res$pearson, 1.0)
  expect_equal(res$spearman, 1.0)
})

test_that("returns correct Pearson for known linear relationship", {
  y_true <- c(1, 2, 3, 4, 5)
  y_pred <- c(2, 4, 6, 8, 10)  # y_pred = 2 * y_true
  res <- evaluate_predictions(y_true, y_pred)
  expect_equal(res$pearson, 1.0)
  expect_equal(res$spearman, 1.0)
})

test_that("negative correlation is captured", {
  y_true <- c(1, 2, 3, 4, 5)
  y_pred <- c(5, 4, 3, 2, 1)
  res <- evaluate_predictions(y_true, y_pred)
  expect_equal(res$pearson, -1.0)
  expect_equal(res$spearman, -1.0)
})

test_that("returns all three metrics", {
  res <- evaluate_predictions(c(1, 2, 3, 4, 5), c(1, 2, 3, 4, 5))
  expect_named(res, c("pearson", "spearman", "ndcg_at_10"))
  expect_true(is.numeric(res$pearson))
  expect_true(is.numeric(res$spearman))
  expect_true(is.numeric(res$ndcg_at_10))
})

test_that("returns NA for single observation", {
  res <- evaluate_predictions(c(5), c(5))
  expect_true(is.na(res$pearson))
  expect_true(is.na(res$spearman))
  expect_true(is.na(res$ndcg_at_10))
})

test_that("lower_is_better trait changes NDCG value for imperfect predictions", {
  # With an imperfect prediction, negating both y_true and y_pred (as
  # lower_is_better does) changes the relevance magnitudes after shifting
  # to positive, producing a different NDCG value.
  y_true <- c(1, 3, 2, 5, 4)
  y_pred <- c(2, 5, 1, 4, 3)  # imperfect prediction

  res_lower <- evaluate_predictions(y_true, y_pred, trait_name = "DTF_blue")
  res_higher <- evaluate_predictions(y_true, y_pred, trait_name = "TGW_blue")

  # Pearson/Spearman are unaffected by the sign flip
  expect_equal(res_lower$pearson, res_higher$pearson)
  expect_equal(res_lower$spearman, res_higher$spearman)

  # NDCG values should differ between the two modes
  expect_false(isTRUE(all.equal(res_lower$ndcg_at_10, res_higher$ndcg_at_10)))

  # Both should still be valid NDCG values in [0, 1]
  expect_true(res_lower$ndcg_at_10 >= 0 && res_lower$ndcg_at_10 <= 1)
  expect_true(res_higher$ndcg_at_10 >= 0 && res_higher$ndcg_at_10 <= 1)
})

test_that("all LOWER_IS_BETTER_TRAITS are recognised", {
  # Use an imperfect prediction so we can verify lower_is_better is active
  # by comparing against a direct calculate_ndcg call
  y_true <- c(1, 3, 2, 5, 4)
  y_pred <- c(2, 5, 1, 4, 3)

  expected <- calculate_ndcg(y_true, y_pred, k = 10, lower_is_better = TRUE)

  for (trait in c("DTF_blue", "DTH_blue", "PtHt_blue")) {
    res <- evaluate_predictions(y_true, y_pred, trait_name = trait)
    expect_equal(res$ndcg_at_10, expected, label = paste("NDCG for", trait))
  }
})

test_that("non-lower_is_better traits don't flip NDCG", {
  # y_pred ranks high values first — correct for higher-is-better -> NDCG = 1
  y_true <- c(5, 4, 3, 2, 1)
  y_pred <- c(5, 4, 3, 2, 1)

  for (trait in c("TGW_blue", "SdLen_blue", "PcleLng_blue", "SdW_z_blue")) {
    res <- evaluate_predictions(y_true, y_pred, trait_name = trait)
    expect_equal(res$ndcg_at_10, 1.0, label = paste("NDCG for", trait))
  }
})

test_that("NULL trait_name defaults to higher-is-better", {
  y_true <- c(5, 4, 3, 2, 1)
  y_pred <- c(5, 4, 3, 2, 1)
  res <- evaluate_predictions(y_true, y_pred, trait_name = NULL)
  expect_equal(res$ndcg_at_10, 1.0)
})
