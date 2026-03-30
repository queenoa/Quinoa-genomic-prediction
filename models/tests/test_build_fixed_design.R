library(testthat)

suppressPackageStartupMessages(
  source(file.path("..", "models", "BayesC", "BayesC_utils.R"))
)

# build_fixed_design requires >= 2 locations (uses model.matrix contrasts).
# This matches the real data which always has multiple locations and years.

test_that("two locations, shared years: location + year-within-loc dummies", {
  # Note: build_fixed_design requires >= 2 year levels globally
  # (line 193 calls model.matrix(~ loc:year) which is dead code but still runs)
  loc  <- c("AUS", "AUS", "PAK", "PAK")
  year <- c("2019", "2020", "2019", "2020")

  Z <- build_fixed_design(loc, year)

  # 1 location dummy (PAK) + 1 year-in-AUS (2020) + 1 year-in-PAK (2020) = 3
  expect_equal(ncol(Z), 3)
  expect_equal(nrow(Z), 4)

  # AUS reference rows have 0 in the location column
  loc_cols <- grep("^loc_", colnames(Z))
  expect_equal(as.numeric(Z[1, loc_cols]), 0)  # AUS 2019
  expect_equal(as.numeric(Z[3, loc_cols]), 1)  # PAK 2019
})

test_that("two locations with different years: location + nested year dummies", {
  loc  <- c("AUS", "AUS", "AUS", "PAK", "PAK")
  year <- c("2017", "2018", "2019", "2019", "2020")

  Z <- build_fixed_design(loc, year)

  # 1 location dummy (PAK) + 2 year-in-AUS (2018, 2019) + 1 year-in-PAK (2020)
  expect_equal(ncol(Z), 4)
  expect_equal(nrow(Z), 5)
})

test_that("dimensions match real AUSPAK structure", {
  # Mimic the actual data: 2 locations, 3+3 years
  loc  <- rep(c("AUS", "AUS", "AUS", "PAK", "PAK", "PAK"), each = 10)
  year <- rep(c("2017", "2018", "2019", "2019", "2020", "2021"), each = 10)

  Z <- build_fixed_design(loc, year)

  # 1 loc dummy + 2 year-in-AUS (2018,2019) + 2 year-in-PAK (2020,2021) = 5
  expect_equal(ncol(Z), 5)
  expect_equal(nrow(Z), 60)
})

test_that("output is a numeric matrix with only 0s and 1s", {
  loc  <- c("AUS", "PAK", "AUS", "PAK")
  year <- c("2019", "2019", "2020", "2020")

  Z <- build_fixed_design(loc, year)

  expect_true(is.matrix(Z))
  expect_true(is.numeric(Z))
  expect_true(all(Z %in% c(0, 1)))
})

test_that("reference location rows have 0 in all location columns", {
  loc  <- c("AUS", "AUS", "PAK", "PAK")
  year <- c("2019", "2020", "2019", "2020")

  Z <- build_fixed_design(loc, year)

  loc_cols <- grep("^loc_", colnames(Z))
  # AUS is the reference level (alphabetically first)
  aus_rows <- which(loc == "AUS")
  expect_true(all(Z[aus_rows, loc_cols] == 0))
})

test_that("each row has at most one 1 in location columns", {
  loc  <- c("AUS", "AUS", "PAK", "PAK")
  year <- c("2017", "2018", "2019", "2020")

  Z <- build_fixed_design(loc, year)

  loc_cols <- grep("^loc_", colnames(Z))
  if (length(loc_cols) > 0) {
    expect_true(all(rowSums(Z[, loc_cols, drop = FALSE]) <= 1))
  }
})

test_that("cbind matrix/data.frame mixing produces correct values with 3+ years", {
  # Regression test: the loop in build_fixed_design starts with a matrix
  # (line 218) then cbinds data.frames (line 221-222). With 3+ years per
  # location the inner loop runs multiple times, exercising the type mix.
  loc  <- c(rep("AUS", 4), rep("PAK", 3))
  year <- c("2017", "2018", "2019", "2020", "2019", "2020", "2021")

  Z <- build_fixed_design(loc, year)

  # Must be a proper numeric matrix, not a data.frame
  expect_true(is.matrix(Z))
  expect_true(is.numeric(Z))
  expect_true(all(Z %in% c(0, 1)))

  # 1 loc (PAK) + 3 year-in-AUS (2018,2019,2020) + 2 year-in-PAK (2020,2021) = 6
  expect_equal(ncol(Z), 6)
  expect_equal(nrow(Z), 7)

  # Verify each year-within-loc dummy is 1 only for the correct row(s)
  # AUS_yr_2018 should be 1 only for row 2 (AUS 2018)
  expect_equal(as.numeric(Z[, "AUS_yr_2018"]), c(0, 1, 0, 0, 0, 0, 0))
  # AUS_yr_2019 should be 1 only for row 3 (AUS 2019)
  expect_equal(as.numeric(Z[, "AUS_yr_2019"]), c(0, 0, 1, 0, 0, 0, 0))
  # AUS_yr_2020 should be 1 only for row 4 (AUS 2020)
  expect_equal(as.numeric(Z[, "AUS_yr_2020"]), c(0, 0, 0, 1, 0, 0, 0))
  # PAK_yr_2020 should be 1 only for row 6 (PAK 2020)
  expect_equal(as.numeric(Z[, "PAK_yr_2020"]), c(0, 0, 0, 0, 0, 1, 0))
  # PAK_yr_2021 should be 1 only for row 7 (PAK 2021)
  expect_equal(as.numeric(Z[, "PAK_yr_2021"]), c(0, 0, 0, 0, 0, 0, 1))
})

test_that("year-within-location dummies are zero outside their location", {
  loc  <- c("AUS", "AUS", "PAK", "PAK")
  year <- c("2017", "2018", "2019", "2020")

  Z <- build_fixed_design(loc, year)

  # Year-in-AUS columns should be 0 for all PAK rows
  aus_yr_cols <- grep("^AUS_yr_", colnames(Z))
  pak_rows <- which(loc == "PAK")
  if (length(aus_yr_cols) > 0 && length(pak_rows) > 0) {
    expect_true(all(Z[pak_rows, aus_yr_cols] == 0))
  }

  # Year-in-PAK columns should be 0 for all AUS rows
  pak_yr_cols <- grep("^PAK_yr_", colnames(Z))
  aus_rows <- which(loc == "AUS")
  if (length(pak_yr_cols) > 0 && length(aus_rows) > 0) {
    expect_true(all(Z[aus_rows, pak_yr_cols] == 0))
  }
})
