library(testthat)
library(Matrix)

suppressPackageStartupMessages(
  source(file.path("..", "GBLUP", "GBLUP_utils.R"))
)

# ── Helper: create mock Ginv sparse matrix with given genotype IDs ──────────

make_mock_ginv <- function(ids) {
  n <- length(ids)
  M <- Matrix(diag(n), sparse = TRUE)
  rownames(M) <- colnames(M) <- ids
  M
}

# ── Helper: create simple phenotype data ────────────────────────────────────

make_pheno <- function(geno_ids, lys = c("AUS_2019", "PAK_2020")) {
  n_ly <- length(lys)
  n_geno <- length(geno_ids)
  data.frame(
    sample.id     = rep(geno_ids, times = n_ly),
    location_year = rep(lys, each = n_geno),
    location      = rep(sub("_.*", "", lys), each = n_geno),
    year          = rep(sub(".*_", "", lys), each = n_geno),
    DTF      = rnorm(n_geno * n_ly),
    stringsAsFactors = FALSE
  )
}

# ============================================================================
# Factor level alignment
# ============================================================================

test_that("sample.id factor levels match G matrix row order", {
  geno_ids <- paste0("G", 1:10)
  Ginv <- make_mock_ginv(geno_ids)
  pheno <- make_pheno(geno_ids[c(3, 5, 7, 1)])  # subset, shuffled

  result <- align_genotypes_to_gmatrix(pheno, Ginv)

  expect_true(is.factor(result$sample.id))
  expect_equal(levels(result$sample.id), geno_ids)
})

test_that("all columns converted to factors", {
  geno_ids <- paste0("G", 1:5)
  Ginv <- make_mock_ginv(geno_ids)
  pheno <- make_pheno(geno_ids)

  result <- align_genotypes_to_gmatrix(pheno, Ginv)

  expect_true(is.factor(result$sample.id))
  expect_true(is.factor(result$location_year))
  expect_true(is.factor(result$location))
  expect_true(is.factor(result$year))
})

test_that("data dimensions are preserved", {
  geno_ids <- paste0("G", 1:8)
  Ginv <- make_mock_ginv(geno_ids)
  pheno <- make_pheno(geno_ids)

  result <- align_genotypes_to_gmatrix(pheno, Ginv)

  expect_equal(nrow(result), nrow(pheno))
  expect_equal(ncol(result), ncol(pheno))
})

test_that("trait values are not modified", {
  set.seed(42)
  geno_ids <- paste0("G", 1:5)
  Ginv <- make_mock_ginv(geno_ids)
  pheno <- make_pheno(geno_ids)

  result <- align_genotypes_to_gmatrix(pheno, Ginv)

  expect_equal(result$DTF, pheno$DTF)
})

# ============================================================================
# G matrix has more genotypes than phenotype data
# ============================================================================

test_that("works when G matrix has extra genotypes", {
  g_ids <- paste0("G", 1:20)          # 20 genotypes in G matrix
  pheno_ids <- paste0("G", c(3, 7, 12))  # only 3 in pheno
  Ginv <- make_mock_ginv(g_ids)
  pheno <- make_pheno(pheno_ids)

  result <- align_genotypes_to_gmatrix(pheno, Ginv)

  # Factor levels should be ALL G matrix genotypes, not just pheno ones
  expect_equal(levels(result$sample.id), g_ids)
  # But data rows should only contain pheno genotypes
  expect_equal(sort(unique(as.character(result$sample.id))),
               sort(pheno_ids))
})

# ============================================================================
# Error: phenotype genotypes not in G matrix
# ============================================================================

test_that("errors when pheno has genotypes not in G matrix", {
  g_ids <- paste0("G", 1:5)
  pheno_ids <- paste0("G", c(1, 2, 99))  # G99 not in G matrix
  Ginv <- make_mock_ginv(g_ids)
  pheno <- make_pheno(pheno_ids)

  expect_error(
    align_genotypes_to_gmatrix(pheno, Ginv),
    "not in G matrix"
  )
})

test_that("error message lists missing genotype IDs", {
  g_ids <- paste0("G", 1:5)
  pheno_ids <- c("G1", "G2", "MISSING1", "MISSING2")
  Ginv <- make_mock_ginv(g_ids)
  pheno <- make_pheno(pheno_ids)

  expect_error(
    align_genotypes_to_gmatrix(pheno, Ginv),
    "MISSING"
  )
})

# ============================================================================
# rowNames attribute fallback (some sparse matrix formats)
# ============================================================================

test_that("reads genotype IDs from rowNames attribute when rownames are NULL", {
  geno_ids <- paste0("G", 1:5)
  Ginv <- make_mock_ginv(geno_ids)

  # Simulate a sparse matrix where rownames() returns NULL but attr has them
  Ginv_no_rownames <- Ginv
  rownames(Ginv_no_rownames) <- NULL
  colnames(Ginv_no_rownames) <- NULL
  attr(Ginv_no_rownames, "rowNames") <- geno_ids

  pheno <- make_pheno(geno_ids)

  result <- align_genotypes_to_gmatrix(pheno, Ginv_no_rownames)
  expect_equal(levels(result$sample.id), geno_ids)
})

# ============================================================================
# Alignment verification (levels(sample.id) must exactly match G rownames)
# ============================================================================

test_that("factor levels are identical to G matrix rownames (same order)", {
  # Reverse order in G matrix
  geno_ids <- paste0("G", 10:1)
  Ginv <- make_mock_ginv(geno_ids)
  pheno <- make_pheno(paste0("G", 1:5))

  result <- align_genotypes_to_gmatrix(pheno, Ginv)

  # Levels should follow G matrix order (10, 9, 8, ..., 1), not natural order
  expect_equal(levels(result$sample.id), geno_ids)
})

# ============================================================================
# Character coercion of sample.id
# ============================================================================

test_that("numeric sample.id is coerced to character then factor", {
  geno_ids <- as.character(101:105)
  Ginv <- make_mock_ginv(geno_ids)

  pheno <- data.frame(
    sample.id     = rep(101:105, 2),  # numeric
    location_year = rep(c("AUS_2019", "PAK_2020"), each = 5),
    location      = rep(c("AUS", "PAK"), each = 5),
    year          = rep(c("2019", "2020"), each = 5),
    DTF      = rnorm(10),
    stringsAsFactors = FALSE
  )

  result <- align_genotypes_to_gmatrix(pheno, Ginv)

  expect_true(is.factor(result$sample.id))
  expect_equal(levels(result$sample.id), geno_ids)
})
