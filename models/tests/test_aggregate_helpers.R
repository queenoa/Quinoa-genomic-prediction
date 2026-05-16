library(testthat)

suppressPackageStartupMessages(
  source(file.path("..", "BayesC", "aggregate_BayesC_results.R"))
)

# ── collect_csvs ─────────────────────────────────────────────────────────────

test_that("collect_csvs returns NULL when no files match", {
  result <- collect_csvs("^NONEXISTENT_PATTERN_.*\\.csv$", "test")
  expect_null(result)
})

test_that("collect_csvs reads and combines matching CSVs", {
  # Create temp CSV files
  tmpdir <- tempdir()
  old_wd <- setwd(tmpdir)
  on.exit(setwd(old_wd))

  df1 <- data.frame(trait = "DTF", pearson = 0.5, stringsAsFactors = FALSE)
  df2 <- data.frame(trait = "DTF", pearson = 0.7, stringsAsFactors = FALSE)
  write.csv(df1, "test_collect_1.csv", row.names = FALSE)
  write.csv(df2, "test_collect_2.csv", row.names = FALSE)

  result <- collect_csvs("^test_collect_[0-9]+\\.csv$", "test")

  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2)
  expect_equal(result$pearson, c(0.5, 0.7))

  # Cleanup
  file.remove("test_collect_1.csv", "test_collect_2.csv")
})

test_that("collect_csvs handles files with different columns via bind_rows", {
  tmpdir <- tempdir()
  old_wd <- setwd(tmpdir)
  on.exit(setwd(old_wd))

  df1 <- data.frame(trait = "DTF", pearson = 0.5, stringsAsFactors = FALSE)
  df2 <- data.frame(trait = "DTF", pearson = 0.7, extra = "x",
                     stringsAsFactors = FALSE)
  write.csv(df1, "test_mixed_1.csv", row.names = FALSE)
  write.csv(df2, "test_mixed_2.csv", row.names = FALSE)

  result <- collect_csvs("^test_mixed_[0-9]+\\.csv$", "test")

  expect_equal(nrow(result), 2)
  # bind_rows fills missing column with NA
  expect_true(is.na(result$extra[1]))
  expect_equal(result$extra[2], "x")

  file.remove("test_mixed_1.csv", "test_mixed_2.csv")
})

# ── compute_summary ──────────────────────────────────────────────────────────

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
    trait = rep("DTF", 6),
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
    trait = rep("DTF", 4),
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
    trait = rep("DTF", 4),
    cv_scheme = rep("CV1", 4),
    location = c("AUS", "AUS", "PAK", "PAK"),
    pearson = c(0.5, 0.6, 0.3, 0.4),
    spearman = c(0.4, 0.5, 0.2, 0.3),
    ndcg_at_10 = c(0.8, 0.9, 0.7, 0.75),
    n_test_genotypes = rep(50, 4),
    stringsAsFactors = FALSE
  )

  result <- compute_summary(cv_data)

  expect_equal(nrow(result), 2)  # AUS and PAK
  aus_row <- result[result$location == "AUS", ]
  pak_row <- result[result$location == "PAK", ]
  expect_equal(aus_row$pearson_mean, 0.55)
  expect_equal(pak_row$pearson_mean, 0.35)
})
