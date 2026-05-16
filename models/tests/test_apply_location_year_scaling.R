library(testthat)

suppressPackageStartupMessages(
  source(file.path("..", "BayesC", "BayesC_utils.R"))
)

# The function requires >= 50 non-NA observations total to apply scaling,
# so test fixtures must meet this threshold.

test_that("z-scores have mean 0 and sd 1 within each location-year", {
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

test_that("NA trait values are preserved after scaling", {
  pheno <- data.frame(
    sample.id     = paste0("G", 1:60),
    location_year = rep("AUS_2019", 60),
    trait1        = c(1:55, rep(NA, 5)),
    stringsAsFactors = FALSE
  )

  scaled <- apply_location_year_scaling(pheno, "trait1")

  # NAs stay NA
  expect_true(all(is.na(scaled$trait1[56:60])))

  # Non-NA values are scaled
  non_na <- scaled$trait1[!is.na(scaled$trait1)]
  expect_equal(mean(non_na), 0, tolerance = 1e-10)
  expect_equal(sd(non_na), 1, tolerance = 1e-10)
})

test_that("zero-variance location-year gets sd forced to 1", {
  # All observed values identical within one location-year -> sd = 0
  # Need >= 50 total non-NA, so use a second location-year with real variance
  pheno <- data.frame(
    sample.id     = paste0("G", 1:100),
    location_year = rep(c("AUS_2019", "PAK_2020"), each = 50),
    trait1        = c(rep(5.0, 50), 1:50),
    stringsAsFactors = FALSE
  )

  scaled <- apply_location_year_scaling(pheno, "trait1")

  # AUS: constant input, sd forced to 1 -> (5 - 5) / 1 = 0 for all
  aus_vals <- scaled$trait1[scaled$location_year == "AUS_2019"]
  expect_true(all(aus_vals == 0))
})

test_that("scaling is independent across location-years", {
  pheno <- data.frame(
    sample.id     = paste0("G", 1:100),
    location_year = rep(c("AUS_2019", "PAK_2020"), each = 50),
    trait1        = c(rep(100, 50), 1:50),  # AUS constant, PAK varies
    stringsAsFactors = FALSE
  )

  scaled <- apply_location_year_scaling(pheno, "trait1")

  # AUS: all constant -> all 0
  aus_vals <- scaled$trait1[scaled$location_year == "AUS_2019"]
  expect_true(all(aus_vals == 0))

  # PAK: should be standard z-scores of 1:50
  pak_vals <- scaled$trait1[scaled$location_year == "PAK_2020"]
  expect_equal(mean(pak_vals), 0, tolerance = 1e-10)
  expect_equal(sd(pak_vals), 1, tolerance = 1e-10)
})

test_that("too few observations triggers a warning and returns unscaled", {
  pheno <- data.frame(
    sample.id     = paste0("G", 1:5),
    location_year = rep("AUS_2019", 5),
    trait1        = c(1:3, NA, NA),  # only 3 non-NA, well below 50
    stringsAsFactors = FALSE
  )

  expect_warning(
    scaled <- apply_location_year_scaling(pheno, "trait1"),
    "Too few observations"
  )

  # Data returned unchanged
  expect_equal(scaled$trait1, pheno$trait1)
})

test_that("ly_mean and ly_sd columns are not left in output", {
  set.seed(42)
  pheno <- data.frame(
    sample.id     = paste0("G", 1:60),
    location_year = rep("AUS_2019", 60),
    trait1        = rnorm(60, 50, 5),
    stringsAsFactors = FALSE
  )

  scaled <- apply_location_year_scaling(pheno, "trait1")
  expect_false("ly_mean" %in% names(scaled))
  expect_false("ly_sd" %in% names(scaled))
})
