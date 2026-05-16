library(testthat)

# Source from RKHS_utils.R (not BayesC_utils.R) to verify the copies are correct
suppressPackageStartupMessages(
  source(file.path("..", "RKHS", "RKHS_utils.R"))
)

# ============================================================================
# VALID_TRAITS and LOWER_IS_BETTER_TRAITS constants
# ============================================================================

test_that("VALID_TRAITS matches BayesC", {
  expected <- c('DTF', 'DTH', 'PtHt', 'PcleLng',
                'SdLen', 'TGW', 'SdW_z')
  expect_equal(VALID_TRAITS, expected)
})

test_that("LOWER_IS_BETTER_TRAITS matches BayesC", {
  expected <- c('DTF', 'DTH', 'PtHt')
  expect_equal(LOWER_IS_BETTER_TRAITS, expected)
})

# ============================================================================
# build_fixed_design (identical to BayesC)
# ============================================================================

test_that("build_fixed_design: two locations, shared years", {
  loc  <- c("AUS", "AUS", "PAK", "PAK")
  year <- c("2019", "2020", "2019", "2020")

  Z <- build_fixed_design(loc, year)

  expect_equal(ncol(Z), 3)
  expect_equal(nrow(Z), 4)
  expect_true(is.matrix(Z))
  expect_true(is.numeric(Z))
  expect_true(all(Z %in% c(0, 1)))
})

test_that("build_fixed_design: year-within-location dummies are zero outside their location", {
  loc  <- c("AUS", "AUS", "PAK", "PAK")
  year <- c("2017", "2018", "2019", "2020")

  Z <- build_fixed_design(loc, year)

  aus_yr_cols <- grep("^AUS_yr_", colnames(Z))
  pak_rows <- which(loc == "PAK")
  if (length(aus_yr_cols) > 0 && length(pak_rows) > 0) {
    expect_true(all(Z[pak_rows, aus_yr_cols] == 0))
  }
})

# ============================================================================
# build_groups
# ============================================================================

test_that("build_groups returns integer factor codes", {
  ly_vec <- c("AUS_2019", "AUS_2020", "PAK_2019", "AUS_2019")
  groups <- build_groups(ly_vec)

  expect_true(is.integer(groups))
  expect_equal(length(groups), 4)
  # Same location-year gets same group
  expect_equal(groups[1], groups[4])
  # Different location-years get different groups
  expect_true(groups[1] != groups[2] || groups[1] != groups[3])
})

test_that("build_groups with single location-year returns all same", {
  groups <- build_groups(rep("AUS_2019", 10))
  expect_true(all(groups == groups[1]))
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

# ============================================================================
# evaluate_predictions
# ============================================================================

test_that("evaluate_predictions returns pearson, spearman, ndcg_at_10", {
  set.seed(42)
  y_true <- rnorm(50)
  y_pred <- y_true + rnorm(50, sd = 0.3)

  res <- evaluate_predictions(y_true, y_pred, trait_name = "DTF")

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

# ============================================================================
# evaluate_per_location_year
# ============================================================================

test_that("evaluate_per_location_year returns list with metrics and predictions", {
  set.seed(42)
  n <- 60
  observed_data <- data.frame(
    sample.id     = paste0("G", 1:n),
    location_year = rep(c("AUS_2019", "PAK_2019", "PAK_2020"), each = 20),
    location      = rep(c("AUS", "PAK", "PAK"), each = 20),
    DTF      = rnorm(n),
    stringsAsFactors = FALSE
  )

  pred_values <- data.frame(
    sample.id     = observed_data$sample.id,
    location_year = observed_data$location_year,
    location      = observed_data$location,
    predicted     = observed_data$DTF + rnorm(n, sd = 0.5),
    stringsAsFactors = FALSE
  )

  result <- evaluate_per_location_year(
    observed_data, pred_values, test_indices = 1:n,
    trait = "DTF", min_genotypes = 5
  )

  expect_true(is.list(result))
  expect_named(result, c("metrics", "predictions"))
  expect_true(is.data.frame(result$metrics))
  expect_true(is.data.frame(result$predictions))
  expect_equal(nrow(result$metrics), 3)  # 3 location-years
  expect_equal(nrow(result$predictions), n)
})

test_that("evaluate_per_location_year: min_genotypes filter works", {
  set.seed(42)
  observed_data <- data.frame(
    sample.id     = paste0("G", 1:5),
    location_year = rep("AUS_2019", 5),
    location      = rep("AUS", 5),
    DTF      = rnorm(5),
    stringsAsFactors = FALSE
  )

  pred_values <- data.frame(
    sample.id     = observed_data$sample.id,
    location_year = observed_data$location_year,
    location      = observed_data$location,
    predicted     = rnorm(5),
    stringsAsFactors = FALSE
  )

  result <- evaluate_per_location_year(
    observed_data, pred_values, test_indices = 1:5,
    trait = "DTF", min_genotypes = 10
  )

  expect_equal(nrow(result$metrics), 0)
})
