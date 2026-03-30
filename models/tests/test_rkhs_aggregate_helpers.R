library(testthat)

suppressPackageStartupMessages(
  source(file.path("..", "models", "RKHS", "aggregate_RKHS_results.R"))
)

# ============================================================================
# read_if_exists
# ============================================================================

test_that("read_if_exists returns NULL for missing file", {
  result <- read_if_exists("nonexistent_file_xyz.csv", "test")
  expect_null(result)
})

test_that("read_if_exists reads existing CSV correctly", {
  tmpfile <- tempfile(fileext = ".csv")
  df <- data.frame(trait = "DTF_blue", pearson = 0.5, stringsAsFactors = FALSE)
  write.csv(df, tmpfile, row.names = FALSE)

  result <- read_if_exists(tmpfile, "test")
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 1)
  expect_equal(result$pearson, 0.5)

  file.remove(tmpfile)
})

test_that("read_if_exists handles corrupt file gracefully", {
  tmpfile <- tempfile(fileext = ".csv")
  writeLines("not,a,valid\ncsv,file,{broken", tmpfile)

  # Should not error, might return the data or NULL depending on how R parses it
  result <- read_if_exists(tmpfile, "test")
  # Main point: no unhandled error
  expect_true(is.null(result) || is.data.frame(result))

  file.remove(tmpfile)
})

# ============================================================================
# compute_summary
# ============================================================================

test_that("compute_summary returns NULL for NULL input", {
  expect_null(compute_summary(NULL))
})

test_that("compute_summary returns NULL for empty data frame", {
  empty_df <- data.frame(
    trait = character(), cv_scheme = character(), location = character(),
    pearson = numeric(), spearman = numeric(), ndcg_at_10 = numeric(),
    n_test_genotypes = integer(), stringsAsFactors = FALSE
  )
  expect_null(compute_summary(empty_df))
})

test_that("compute_summary groups by trait, cv_scheme, location", {
  cv_data <- data.frame(
    trait = rep("DTF_blue", 6),
    cv_scheme = rep(c("CV1", "CV2"), each = 3),
    location = rep("AUS", 6),
    pearson = c(0.5, 0.6, 0.7, 0.3, 0.4, 0.5),
    spearman = c(0.4, 0.5, 0.6, 0.2, 0.3, 0.4),
    ndcg_at_10 = c(0.8, 0.85, 0.9, 0.7, 0.75, 0.8),
    n_test_genotypes = rep(50, 6),
    stringsAsFactors = FALSE
  )

  result <- compute_summary(cv_data)

  expect_equal(nrow(result), 2)  # CV1 and CV2
  expect_true(all(c("pearson_mean", "pearson_sd", "pearson_min", "pearson_max",
                     "spearman_mean", "spearman_sd",
                     "ndcg_at_10_mean", "ndcg_at_10_sd",
                     "n_evaluations", "mean_n_test_genotypes")
                   %in% names(result)))
})

test_that("compute_summary calculates correct mean and sd", {
  cv_data <- data.frame(
    trait = rep("DTF_blue", 4),
    cv_scheme = rep("CV1", 4),
    location = rep("AUS", 4),
    pearson = c(0.4, 0.6, 0.8, 1.0),
    spearman = c(0.3, 0.5, 0.7, 0.9),
    ndcg_at_10 = c(0.7, 0.8, 0.9, 1.0),
    n_test_genotypes = rep(50, 4),
    stringsAsFactors = FALSE
  )

  result <- compute_summary(cv_data)

  expect_equal(result$pearson_mean, mean(c(0.4, 0.6, 0.8, 1.0)))
  expect_equal(result$pearson_sd, sd(c(0.4, 0.6, 0.8, 1.0)))
  expect_equal(result$pearson_min, 0.4)
  expect_equal(result$pearson_max, 1.0)
  expect_equal(result$n_evaluations, 4)
})

test_that("compute_summary separates locations", {
  cv_data <- data.frame(
    trait = rep("DTF_blue", 4),
    cv_scheme = rep("CV1", 4),
    location = c("AUS", "AUS", "PAK", "PAK"),
    pearson = c(0.5, 0.6, 0.3, 0.4),
    spearman = c(0.4, 0.5, 0.2, 0.3),
    ndcg_at_10 = c(0.8, 0.9, 0.7, 0.75),
    n_test_genotypes = rep(50, 4),
    stringsAsFactors = FALSE
  )

  result <- compute_summary(cv_data)

  expect_equal(nrow(result), 2)
  aus_row <- result[result$location == "AUS", ]
  pak_row <- result[result$location == "PAK", ]
  expect_equal(aus_row$pearson_mean, 0.55)
  expect_equal(pak_row$pearson_mean, 0.35)
})

# ============================================================================
# CrossLoc scheme names are preserved correctly
# ============================================================================

test_that("compute_summary handles CrossLoc scheme names", {
  cv_data <- data.frame(
    trait = rep("DTF_blue", 4),
    cv_scheme = c("CrossLoc_AUS->PAK", "CrossLoc_AUS->PAK",
                  "CrossLoc_PAK->AUS", "CrossLoc_PAK->AUS"),
    location = c("PAK", "PAK", "AUS", "AUS"),
    pearson = c(0.3, 0.4, 0.5, 0.6),
    spearman = c(0.2, 0.3, 0.4, 0.5),
    ndcg_at_10 = c(0.6, 0.7, 0.8, 0.9),
    n_test_genotypes = rep(50, 4),
    stringsAsFactors = FALSE
  )

  result <- compute_summary(cv_data)

  expect_equal(nrow(result), 2)
  expect_true(all(grepl("CrossLoc_", result$cv_scheme)))
})
