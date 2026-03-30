library(testthat)

suppressPackageStartupMessages(
  source(file.path("..", "models", "BayesC", "BayesC_utils.R"))
)

# ── build_obs_marker_matrix ──────────────────────────────────────────────────

test_that("subsets marker matrix to requested sample IDs in order", {
  X <- matrix(1:12, nrow = 4, ncol = 3,
              dimnames = list(c("G1", "G2", "G3", "G4"), c("S1", "S2", "S3")))

  result <- build_obs_marker_matrix(c("G2", "G4"), X)

  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 3)
  expect_equal(as.numeric(result[1, ]), as.numeric(X["G2", ]))
  expect_equal(as.numeric(result[2, ]), as.numeric(X["G4", ]))
})

test_that("repeated sample IDs produce repeated rows", {
  X <- matrix(1:6, nrow = 2, ncol = 3,
              dimnames = list(c("G1", "G2"), c("S1", "S2", "S3")))

  # G1 appears in multiple environments -> repeated in observation-level matrix
  result <- build_obs_marker_matrix(c("G1", "G2", "G1"), X)

  expect_equal(nrow(result), 3)
  expect_equal(as.numeric(result[1, ]), as.numeric(X["G1", ]))
  expect_equal(as.numeric(result[3, ]), as.numeric(X["G1", ]))
})

test_that("rownames are removed from output", {
  X <- matrix(1:6, nrow = 2, ncol = 3,
              dimnames = list(c("G1", "G2"), c("S1", "S2", "S3")))

  result <- build_obs_marker_matrix(c("G1", "G2"), X)
  expect_null(rownames(result))
})

test_that("column names are preserved", {
  X <- matrix(1:6, nrow = 2, ncol = 3,
              dimnames = list(c("G1", "G2"), c("SNP1", "SNP2", "SNP3")))

  result <- build_obs_marker_matrix(c("G1"), X)
  expect_equal(colnames(result), c("SNP1", "SNP2", "SNP3"))
})

# ── build_groups ─────────────────────────────────────────────────────────────

test_that("returns consecutive integers starting from 1", {
  ly_vec <- c("AUS_2017", "AUS_2017", "PAK_2019", "PAK_2019", "AUS_2018")
  groups <- build_groups(ly_vec)

  expect_true(is.integer(groups))
  expect_equal(length(groups), 5)
  expect_equal(min(groups), 1L)
  expect_equal(max(groups), length(unique(ly_vec)))
})

test_that("same location-year gets the same group integer", {
  ly_vec <- c("AUS_2017", "PAK_2019", "AUS_2017", "PAK_2019")
  groups <- build_groups(ly_vec)

  expect_equal(groups[1], groups[3])  # both AUS_2017
  expect_equal(groups[2], groups[4])  # both PAK_2019
  expect_false(groups[1] == groups[2])  # different location-years
})

test_that("single location-year returns all 1s", {
  groups <- build_groups(rep("AUS_2019", 5))
  expect_true(all(groups == 1L))
})
