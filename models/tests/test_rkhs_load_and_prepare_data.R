library(testthat)

suppressPackageStartupMessages(
  source(file.path("..", "models", "RKHS", "RKHS_utils.R"))
)

MARKER_FILE <- file.path("..", "data", "AUSPAK_test_subset_1k.raw")
PHENO_FILE  <- file.path("..", "data", "AUSPAK_phenotypes_means_BLUEs.csv")

# ── Input validation ─────────────────────────────────────────────────────────

test_that("invalid trait name is rejected", {
  expect_error(
    load_and_prepare_data("NOT_A_TRAIT", MARKER_FILE, PHENO_FILE),
    "Invalid trait"
  )
})

test_that("missing phenotype file is rejected", {
  expect_error(
    load_and_prepare_data("DTF_blue", MARKER_FILE, "nonexistent.csv"),
    "Phenotype file not found"
  )
})

# ── Successful load with test data ───────────────────────────────────────────

test_that("returns a list with pheno, K_geno_list, h_values, med_D", {
  # Use a temp directory for kernel checkpoint to avoid clobbering real one
  tmpdir <- tempdir()
  kcp <- file.path(tmpdir, "test_kernels.RData")
  if (file.exists(kcp)) file.remove(kcp)

  dat <- load_and_prepare_data("DTF_blue", MARKER_FILE, PHENO_FILE,
                                kernel_checkpoint = kcp)

  expect_true(is.list(dat))
  expect_named(dat, c("pheno", "K_geno_list", "h_values", "med_D"))
  expect_true(is.data.frame(dat$pheno))
  expect_true(is.list(dat$K_geno_list))
  expect_equal(length(dat$K_geno_list), 3)
  expect_equal(length(dat$h_values), 3)
  expect_true(is.numeric(dat$med_D))

  # Cleanup
  if (file.exists(kcp)) file.remove(kcp)
})

test_that("genotype-level kernels are square and match genotype count", {
  tmpdir <- tempdir()
  kcp <- file.path(tmpdir, "test_kernels2.RData")
  if (file.exists(kcp)) file.remove(kcp)

  dat <- load_and_prepare_data("DTF_blue", MARKER_FILE, PHENO_FILE,
                                kernel_checkpoint = kcp)

  n_geno <- length(unique(dat$pheno$sample.id))
  for (i in 1:3) {
    K <- dat$K_geno_list[[i]]
    expect_true(is.matrix(K))
    expect_equal(nrow(K), ncol(K))
    # Kernel dimensions should equal number of genotypes in marker file
    # (not necessarily same as pheno genotypes if some markers are missing)
    expect_gte(nrow(K), n_geno)
  }

  if (file.exists(kcp)) file.remove(kcp)
})

# ── Kernel checkpoint ────────────────────────────────────────────────────────

test_that("kernel checkpoint is created and loadable", {
  tmpdir <- tempdir()
  kcp <- file.path(tmpdir, "test_kernels_ckpt.RData")
  if (file.exists(kcp)) file.remove(kcp)

  dat <- load_and_prepare_data("DTF_blue", MARKER_FILE, PHENO_FILE,
                                kernel_checkpoint = kcp)

  expect_true(file.exists(kcp))

  # Verify checkpoint contents
  env <- new.env()
  load(kcp, envir = env)
  expect_true(exists("K_geno_list", envir = env))
  expect_true(exists("geno_ids_kernel", envir = env))
  expect_true(exists("h_values", envir = env))
  expect_true(exists("med_D", envir = env))

  if (file.exists(kcp)) file.remove(kcp)
})

test_that("loading from checkpoint produces identical results", {
  tmpdir <- tempdir()
  kcp <- file.path(tmpdir, "test_kernels_reload.RData")
  if (file.exists(kcp)) file.remove(kcp)

  # First call: compute kernels
  dat1 <- load_and_prepare_data("DTF_blue", MARKER_FILE, PHENO_FILE,
                                 kernel_checkpoint = kcp)
  expect_true(file.exists(kcp))

  # Second call: load from checkpoint
  dat2 <- load_and_prepare_data("DTF_blue", MARKER_FILE, PHENO_FILE,
                                 kernel_checkpoint = kcp)

  # Kernels should be identical
  for (i in 1:3) {
    expect_equal(dat1$K_geno_list[[i]], dat2$K_geno_list[[i]],
                 tolerance = 1e-15,
                 info = paste("Kernel", i, "should be identical from checkpoint"))
  }
  expect_equal(dat1$h_values, dat2$h_values)
  expect_equal(dat1$med_D, dat2$med_D)

  if (file.exists(kcp)) file.remove(kcp)
})

# ── Phenotype alignment and preprocessing ────────────────────────────────────

test_that("all pheno sample IDs exist in kernel matrices", {
  tmpdir <- tempdir()
  kcp <- file.path(tmpdir, "test_kernels_align.RData")
  if (file.exists(kcp)) file.remove(kcp)

  dat <- load_and_prepare_data("DTF_blue", MARKER_FILE, PHENO_FILE,
                                kernel_checkpoint = kcp)

  pheno_ids <- unique(dat$pheno$sample.id)
  kernel_ids <- rownames(dat$K_geno_list[[1]])
  expect_true(all(pheno_ids %in% kernel_ids))

  if (file.exists(kcp)) file.remove(kcp)
})

test_that("pheno contains required columns", {
  tmpdir <- tempdir()
  kcp <- file.path(tmpdir, "test_kernels_cols.RData")
  if (file.exists(kcp)) file.remove(kcp)

  dat <- load_and_prepare_data("DTF_blue", MARKER_FILE, PHENO_FILE,
                                kernel_checkpoint = kcp)
  required <- c("sample.id", "location", "year", "location_year", "DTF_blue")
  expect_true(all(required %in% names(dat$pheno)))

  if (file.exists(kcp)) file.remove(kcp)
})

test_that("trait values are z-scored (mean ~0 per location-year)", {
  tmpdir <- tempdir()
  kcp <- file.path(tmpdir, "test_kernels_zscore.RData")
  if (file.exists(kcp)) file.remove(kcp)

  dat <- load_and_prepare_data("DTF_blue", MARKER_FILE, PHENO_FILE,
                                kernel_checkpoint = kcp)
  pheno <- dat$pheno

  for (ly in unique(pheno$location_year)) {
    vals <- pheno$DTF_blue[pheno$location_year == ly & !is.na(pheno$DTF_blue)]
    if (length(vals) > 10) {
      expect_equal(mean(vals), 0, tolerance = 1e-8,
                   label = paste(ly, "mean"))
    }
  }

  if (file.exists(kcp)) file.remove(kcp)
})

# ── Empty location-years NOT dropped (RKHS design decision) ─────────────────

test_that("empty location-years are NOT dropped (intentional for RKHS)", {
  tmpdir <- tempdir()
  kcp <- file.path(tmpdir, "test_kernels_empty.RData")
  if (file.exists(kcp)) file.remove(kcp)

  # SdW_z_blue often has empty location-years
  dat <- load_and_prepare_data("SdW_z_blue", MARKER_FILE, PHENO_FILE,
                                kernel_checkpoint = kcp)
  pheno <- dat$pheno

  # Check if any location-year has 0 observed values
  ly_obs_counts <- tapply(!is.na(pheno$SdW_z_blue), pheno$location_year, sum)
  has_empty <- any(ly_obs_counts == 0)

  # Either there are no empty LYs in this trait, or they are kept
  # (Unlike BayesC which drops them)
  # The key test: total rows should include all location-years from phenotype file
  original_pheno <- read.csv(PHENO_FILE, stringsAsFactors = FALSE)
  original_pheno <- original_pheno[original_pheno$sample.id %in%
                                     rownames(dat$K_geno_list[[1]]), ]
  original_lys <- unique(original_pheno$location_year)
  loaded_lys <- unique(pheno$location_year)

  # RKHS keeps all location-years (including potentially empty ones)
  expect_true(all(original_lys %in% loaded_lys),
              info = "RKHS should keep all location-years, even empty ones")

  if (file.exists(kcp)) file.remove(kcp)
})

# ── All traits loadable ─────────────────────────────────────────────────────

test_that("all valid traits can be loaded", {
  tmpdir <- tempdir()
  kcp <- file.path(tmpdir, "test_kernels_alltraits.RData")
  if (file.exists(kcp)) file.remove(kcp)

  for (trait in VALID_TRAITS) {
    dat <- load_and_prepare_data(trait, MARKER_FILE, PHENO_FILE,
                                  kernel_checkpoint = kcp)
    expect_true(is.data.frame(dat$pheno), label = paste("pheno for", trait))
    expect_equal(length(dat$K_geno_list), 3, label = paste("kernels for", trait))
  }

  if (file.exists(kcp)) file.remove(kcp)
})
