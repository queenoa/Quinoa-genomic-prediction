library(testthat)

suppressPackageStartupMessages(
  source(file.path("..", "BayesC", "BayesC_utils.R"))
)

MARKER_FILE <- file.path("..", "..", "data", "AUSPAK_test_subset_1k.raw")
PHENO_FILE  <- file.path("..", "..", "data", "AUSPAK_phenotypes_GP_input.csv")

# ── Input validation ─────────────────────────────────────────────────────────

test_that("invalid trait name is rejected", {
  expect_error(
    load_and_prepare_data("NOT_A_TRAIT", MARKER_FILE, PHENO_FILE),
    "Invalid trait"
  )
})

test_that("missing marker file is rejected", {
  expect_error(
    load_and_prepare_data("DTF", "nonexistent.raw", PHENO_FILE),
    "Marker file not found"
  )
})

test_that("missing phenotype file is rejected", {
  expect_error(
    load_and_prepare_data("DTF", MARKER_FILE, "nonexistent.csv"),
    "Phenotype file not found"
  )
})

# ── Successful load with test data ───────────────────────────────────────────

test_that("returns a list with pheno and X_geno", {
  dat <- load_and_prepare_data("DTF", MARKER_FILE, PHENO_FILE)

  expect_true(is.list(dat))
  expect_named(dat, c("pheno", "X_geno"))
  expect_true(is.data.frame(dat$pheno))
  expect_true(is.matrix(dat$X_geno))
})

test_that("marker matrix has no NAs after imputation", {
  dat <- load_and_prepare_data("DTF", MARKER_FILE, PHENO_FILE)
  expect_equal(sum(is.na(dat$X_geno)), 0)
})

test_that("marker matrix has 1000 SNP columns (test subset)", {
  dat <- load_and_prepare_data("DTF", MARKER_FILE, PHENO_FILE)
  expect_equal(ncol(dat$X_geno), 1000)
})

test_that("all phenotype sample IDs exist in marker matrix", {
  dat <- load_and_prepare_data("DTF", MARKER_FILE, PHENO_FILE)
  expect_true(all(dat$pheno$sample.id %in% rownames(dat$X_geno)))
})

test_that("pheno sample IDs can be used to index X_geno without reordering", {
  # Critical: build_obs_marker_matrix uses pheno$sample.id to subset X_geno.
  # If a sample.id appears in pheno but not in X_geno rownames, the subsetting
  # would fail or return wrong rows silently.
  dat <- load_and_prepare_data("DTF", MARKER_FILE, PHENO_FILE)
  ids <- dat$pheno$sample.id

  # Every ID in pheno must be a valid rowname in X_geno
  expect_true(all(ids %in% rownames(dat$X_geno)))

  # Subsetting should return the correct number of rows (including repeats
  # for genotypes appearing in multiple location-years)
  X_sub <- dat$X_geno[ids, , drop = FALSE]
  expect_equal(nrow(X_sub), nrow(dat$pheno))
})

test_that("pheno contains required columns", {
  dat <- load_and_prepare_data("DTF", MARKER_FILE, PHENO_FILE)
  required <- c("sample.id", "location", "year", "location_year", "DTF")
  expect_true(all(required %in% names(dat$pheno)))
})

test_that("trait values are z-scored (approximately mean 0 per location-year)", {
  dat <- load_and_prepare_data("DTF", MARKER_FILE, PHENO_FILE)
  pheno <- dat$pheno

  for (ly in unique(pheno$location_year)) {
    vals <- pheno$DTF[pheno$location_year == ly & !is.na(pheno$DTF)]
    if (length(vals) > 10) {
      expect_equal(mean(vals), 0, tolerance = 1e-8,
                   label = paste(ly, "mean"))
    }
  }
})

test_that("location-years with zero observed values are dropped", {
  # Verify that any location-year with 0 observed values is removed
  dat <- load_and_prepare_data("SdW_z", MARKER_FILE, PHENO_FILE)
  pheno <- dat$pheno

  for (ly in unique(pheno$location_year)) {
    n_obs <- sum(!is.na(pheno$SdW_z[pheno$location_year == ly]))
    expect_gt(n_obs, 0, label = paste(ly, "should have > 0 observed values"))
  }
})

test_that("all valid traits can be loaded", {
  for (trait in VALID_TRAITS) {
    dat <- load_and_prepare_data(trait, MARKER_FILE, PHENO_FILE)
    expect_true(is.data.frame(dat$pheno), label = paste("pheno for", trait))
    expect_true(is.matrix(dat$X_geno), label = paste("X_geno for", trait))
  }
})
