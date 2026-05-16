library(testthat)

# Source the utils
suppressPackageStartupMessages(
  source(file.path("..", "BayesC", "BayesC_utils.R"))
)

test_that("perfect ranking returns NDCG = 1", {
  y_true <- c(10, 8, 6, 4, 2)
  y_pred <- c(10, 8, 6, 4, 2)  # same order
  expect_equal(calculate_ndcg(y_true, y_pred, k = 5), 1.0)
})

test_that("perfect ranking with k < n returns NDCG = 1", {
  y_true <- c(10, 8, 6, 4, 2)
  y_pred <- c(10, 8, 6, 4, 2)
  expect_equal(calculate_ndcg(y_true, y_pred, k = 3), 1.0)
})

test_that("worst (reversed) ranking returns NDCG < 1", {
  y_true <- c(10, 8, 6, 4, 2)
  y_pred <- c(2, 4, 6, 8, 10)  # reversed
  result <- calculate_ndcg(y_true, y_pred, k = 5)
  expect_lt(result, 1.0)
  expect_gt(result, 0.0)
})

test_that("k is clamped to length of input", {
  y_true <- c(5, 3, 1)
  y_pred <- c(5, 3, 1)
  # k=10 but only 3 elements — should still work and return 1.0
  expect_equal(calculate_ndcg(y_true, y_pred, k = 10), 1.0)
})

test_that("single element returns 1.0", {
  expect_equal(calculate_ndcg(c(5), c(5), k = 10), 1.0)
})

test_that("lower_is_better flips the ranking", {
  # For lower-is-better: true best = smallest value
  y_true <- c(1, 2, 3, 4, 5)
  y_pred <- c(1, 2, 3, 4, 5)  # pred correctly ranks 1 as "best" (lowest)

  # With lower_is_better=TRUE, the predicted ranking of the best (lowest)
  # items should give NDCG = 1
  result <- calculate_ndcg(y_true, y_pred, k = 5, lower_is_better = TRUE)
  expect_equal(result, 1.0)
})

test_that("lower_is_better penalises ranking that favours high values", {
  y_true <- c(1, 2, 3, 4, 5)
  # Prediction ranks high values as best — wrong for lower-is-better
  y_pred <- c(5, 4, 3, 2, 1)
  result <- calculate_ndcg(y_true, y_pred, k = 5, lower_is_better = TRUE)
  expect_lt(result, 1.0)
})

test_that("tied true values give NDCG = 1 regardless of pred order", {
  # All true values identical — every ranking is equally "correct",
  # so NDCG should be 1.0 no matter what y_pred looks like
  y_true <- c(5, 5, 5, 5)

  expect_equal(calculate_ndcg(y_true, c(1, 2, 3, 4), k = 4), 1.0)
  expect_equal(calculate_ndcg(y_true, c(4, 3, 2, 1), k = 4), 1.0)
  expect_equal(calculate_ndcg(y_true, c(9, 1, 9, 1), k = 4), 1.0)
})

test_that("NDCG handles negative values via shifting", {
  y_true <- c(-2, -1, 0, 1, 2)
  y_pred <- c(-2, -1, 0, 1, 2)
  result <- calculate_ndcg(y_true, y_pred, k = 5)
  expect_equal(result, 1.0)
})

test_that("known NDCG@3 value for a specific ranking", {
  # True relevance: [10, 8, 6, 4, 2]
  # Predicted order puts items 3rd, 1st, 2nd, 5th, 4th best on top
  y_true <- c(10, 8, 6, 4, 2)
  y_pred <- c(3,  9, 7, 1, 5)
  # Predicted top-3 (by y_pred desc): indices 2,3,5 -> true vals 8,6,2
  # DCG@3 = 8/log2(2) + 6/log2(3) + 2/log2(4)
  # IDCG@3 = 10/log2(2) + 8/log2(3) + 6/log2(4)
  expected_dcg  <- 8/log2(2) + 6/log2(3) + 2/log2(4)
  expected_idcg <- 10/log2(2) + 8/log2(3) + 6/log2(4)
  expected_ndcg <- expected_dcg / expected_idcg
  result <- calculate_ndcg(y_true, y_pred, k = 3)
  expect_equal(result, expected_ndcg, tolerance = 1e-6)
})
