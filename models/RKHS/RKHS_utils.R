################################################################################
# RKHS_utils.R — Shared utility functions for RKHS multi-environment CV
#
# Sourced by:
#   RKHS_CV1_single_iter.R
#   RKHS_CV2_single_iter.R
#   RKHS_CV0_CrossLoc.R
#
# Contains:
#   - Data loading, preprocessing, and kernel computation
#   - Z-score standardisation (identical to BayesC/GBLUP pipelines)
#   - Observation-level kernel expansion
#   - Fixed-effect design matrix (location + year-within-location)
#   - BGLR RKHS model fitting wrapper (multi-kernel)
#   - Evaluation metrics (Pearson, Spearman, NDCG@10)
#   - Per-location-year evaluation helper
#
# Method:
#   Gaussian kernel with multi-kernel averaging (KA):
#   → Markers are centered and scaled before distance computation
#   → K(xi, xj) = exp(-h * D_ij)  where D is squared Euclidean distance
#   → D is normalised by mean(D)
#   → Three kernels at h = 1/median(D) * {1/5, 1, 5}, following BGLR paper
#     (Pérez & de los Campos 2014, Box 11)
#   → BGLR estimates the variance component for each kernel, effectively
#     performing Bayesian model averaging over bandwidths
#
# Reference:
#   Pérez & de los Campos (2014) Genetics 198:483-495
#   de los Campos et al. (2010) Genetics Research 92:295-308
#   Cuevas et al. (2016) Plant Genome 9(3)
################################################################################

library(BGLR)
library(dplyr)
library(data.table)

# ── Global constants ─────────────────────────────────────────────────────────

VALID_TRAITS <- c('DTF', 'DTH', 'PtHt', 'PcleLng',
                  'SdLen', 'TGW', 'SdW_z')

LOWER_IS_BETTER_TRAITS <- c('DTF', 'DTH', 'PtHt')

# ============================================================================
# DATA LOADING, KERNEL COMPUTATION, AND PREPROCESSING
# ============================================================================

#' Load phenotype and kernel data for RKHS multi-environment CV.
#'
#' The genotype-level kernels (n_geno × n_geno) are computed once from the
#' full marker matrix and checkpointed to RKHS_kernels.RData. On subsequent
#' calls the checkpoint is loaded, avoiding redundant computation.
#'
#' Unlike BayesC, empty location-years (100% missing for the target trait) are
#' NOT dropped. The RKHS kernel is precomputed at the genotype level and
#' expanded to observation level by simple indexing — the computational cost
#' of extra NA rows is negligible, and keeping them matches GBLUP behaviour.
#'
#' @param trait Character: trait name (validated against VALID_TRAITS)
#' @param marker_file Character: path to plink --export A .raw file
#' @param pheno_file Character: path to phenotype CSV
#' @param kernel_checkpoint Character: path for kernel checkpoint RData file
#' @return List with components: pheno, K_geno_list, h_values, med_D

load_and_prepare_data <- function(trait, marker_file,
                                   pheno_file = "../AUSPAK_phenotypes_GP_input.csv",
                                   kernel_checkpoint = "RKHS_kernels.RData") {

  # ── Validate inputs ──────────────────────────────────────────────────────
  if (!trait %in% VALID_TRAITS) {
    stop("Invalid trait: '", trait, "'\n",
         "  Valid traits: ", paste(VALID_TRAITS, collapse = ", "))
  }
  if (!file.exists(pheno_file)) {
    stop("Phenotype file not found: ", pheno_file)
  }

  # ── Load phenotype data ──────────────────────────────────────────────────
  pheno <- read.csv(pheno_file, stringsAsFactors = FALSE)

  required_cols <- c("sample.id", "location", "year", "location_year", trait)
  missing_cols <- setdiff(required_cols, names(pheno))
  if (length(missing_cols) > 0) {
    stop("Missing columns in phenotype file: ", paste(missing_cols, collapse = ", "))
  }

  cat("Phenotype data:", nrow(pheno), "rows,", ncol(pheno), "columns\n")
  cat("Location-years:", paste(sort(unique(pheno$location_year)), collapse = ", "), "\n")

  # ── Load or compute genotype-level kernels ───────────────────────────────
  if (file.exists(kernel_checkpoint)) {

    cat("Loading pre-computed kernels from", kernel_checkpoint, "...\n")
    load(kernel_checkpoint)
    cat("Loaded:", length(geno_ids_kernel), "genotypes,",
        length(K_geno_list), "kernels (h =",
        paste(round(h_values, 5), collapse = ", "), ")\n")

  } else {

    if (!file.exists(marker_file)) {
      stop("Marker file not found: ", marker_file)
    }

    cat("Computing kernels from", marker_file, "...\n")

    raw <- fread(marker_file, header = TRUE, check.names = FALSE)
    geno_ids_kernel <- as.character(raw$IID)
    X <- as.matrix(raw[, -(1:6)])
    rownames(X) <- geno_ids_kernel
    rm(raw); gc()

    cat("Marker matrix:", nrow(X), "genotypes x", ncol(X), "SNPs\n")

    # Mean-impute missing genotypes before scaling
    na_count <- sum(is.na(X))
    cat("Missing genotypes:", na_count, "(",
        round(mean(is.na(X)) * 100, 4), "% of all marker calls)\n")

    if (na_count > 0) {
      imp_vals <- colMeans(X, na.rm = TRUE)
      ix <- which(is.na(X), arr.ind = TRUE)
      X[ix] <- imp_vals[ix[, 2]]
      cat("Mean-imputed", nrow(ix), "missing calls across",
          length(unique(ix[, 2])), "markers.\n")
    }
    stopifnot(sum(is.na(X)) == 0)

    # Center and scale so each SNP contributes equally regardless of MAF
    cat("Centering and scaling marker matrix...\n")
    X <- scale(X, center = TRUE, scale = TRUE)

    # Remove monomorphic markers (zero variance → NaN after scaling)
    nan_cols <- which(colSums(is.nan(X)) > 0)
    if (length(nan_cols) > 0) {
      cat("Removing", length(nan_cols), "monomorphic markers after scaling.\n")
      X <- X[, -nan_cols]
    }
    cat("Markers after QC:", ncol(X), "\n")

    # Compute squared Euclidean distance, normalise by mean
    cat("Computing squared Euclidean distance matrix...\n")
    D <- as.matrix(dist(X, method = "euclidean"))^2
    D <- D / mean(D)
    cat("Distance matrix: normalised by mean(D)\n")

    # Bandwidth: base h = 1/median(D), then three kernels spanning 25-fold range
    med_D <- median(D[lower.tri(D)])
    h_base <- 1 / med_D
    h_values <- h_base * c(1/5, 1, 5)

    cat("Median distance (off-diag):", round(med_D, 4), "\n")
    cat("Bandwidth values (h):", round(h_values, 5), "\n")

    # Pre-compute the three genotype-level kernels
    K_geno_list <- lapply(h_values, function(h) {
      K <- exp(-h * D)
      rownames(K) <- colnames(K) <- geno_ids_kernel
      K
    })
    names(K_geno_list) <- paste0("h_", round(h_values, 5))

    cat("Kernels computed. Dimensions:", nrow(K_geno_list[[1]]), "x",
        ncol(K_geno_list[[1]]), "\n")

    # Free memory
    rm(X, D); gc()

    # Save checkpoint
    save(K_geno_list, geno_ids_kernel, h_values, med_D,
         file = kernel_checkpoint)
    cat("Kernel checkpoint saved to", kernel_checkpoint, "\n")
  }

  # ── Align phenotype and genotype data ────────────────────────────────────
  pheno$sample.id     <- as.character(pheno$sample.id)
  pheno$location      <- as.character(pheno$location)
  pheno$year          <- as.character(pheno$year)
  pheno$location_year <- as.character(pheno$location_year)

  pheno <- pheno[pheno$sample.id %in% geno_ids_kernel, ]
  cat("Observations with genotype data:", nrow(pheno), "\n")
  cat("Unique genotypes:", length(unique(pheno$sample.id)), "\n")

  stopifnot(all(pheno$sample.id %in% rownames(K_geno_list[[1]])))

  # Z-score standardisation by location-year
  cat("Applying z-score standardisation by location-year...\n")
  pheno <- apply_location_year_scaling(pheno, trait)

  # Report missingness
  n_obs <- sum(!is.na(pheno[[trait]]))
  n_total <- nrow(pheno)
  cat("Trait", trait, ":", n_obs, "observed values out of", n_total,
      "(", round(100 * (n_total - n_obs) / n_total, 1), "% missing)\n")

  # Per location-year breakdown
  cat("\nPer location-year breakdown:\n")
  ly_summary <- pheno %>%
    group_by(location_year) %>%
    summarise(n_total = n(),
              n_observed = sum(!is.na(.data[[trait]])),
              pct_missing = round(100 * mean(is.na(.data[[trait]])), 1),
              .groups = 'drop')
  print(as.data.frame(ly_summary))
  cat("\n")

  # NOTE: Unlike BayesC, we do NOT drop empty location-years.
  # Kernel expansion is a cheap indexing operation, so extra NA rows
  # have negligible cost. This matches GBLUP behaviour.

  return(list(
    pheno        = pheno,
    K_geno_list  = K_geno_list,
    h_values     = h_values,
    med_D        = med_D
  ))
}


# ============================================================================
# Z-SCORE STANDARDISATION (identical to BayesC_utils.R)
# ============================================================================

apply_location_year_scaling <- function(pheno_data, trait) {
  scaled_data <- pheno_data
  trait_data <- pheno_data[!is.na(pheno_data[[trait]]), ]

  if (nrow(trait_data) < 50) {
    warning(paste("Too few observations for z-score scaling of", trait))
    return(scaled_data)
  }

  ly_stats <- trait_data %>%
    group_by(location_year) %>%
    summarise(
      ly_mean = mean(.data[[trait]], na.rm = TRUE),
      ly_sd   = sd(.data[[trait]], na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(ly_sd = ifelse(ly_sd == 0 | is.na(ly_sd), 1, ly_sd))

  scaled_data <- scaled_data %>%
    left_join(ly_stats, by = "location_year") %>%
    mutate(
      !!trait := ifelse(!is.na(.data[[trait]]),
                        (.data[[trait]] - ly_mean) / ly_sd,
                        NA)
    ) %>%
    select(-ly_mean, -ly_sd)

  return(scaled_data)
}


# ============================================================================
# OBSERVATION-LEVEL CONSTRUCTION HELPERS
# ============================================================================

#' Expand genotype-level kernel to observation-level kernel.
#'
#' For multi-environment data in long format, observation i corresponds to
#' genotype g_i. The observation-level kernel is:
#'   K_obs[i,j] = K_geno[g_i, g_j]
#'
#' This is the kernel analogue of build_obs_marker_matrix() in BayesC,
#' but much cheaper — just indexing, no matrix multiplication.
#'
#' @param sample_ids Character vector of sample IDs (length = n_obs),
#'        in the order they appear in the model data.
#' @param K_geno Square kernel matrix with genotype IDs as row/colnames.
#' @return Square observation-level kernel matrix (n_obs × n_obs).

build_obs_kernel <- function(sample_ids, K_geno) {
  K_obs <- K_geno[sample_ids, sample_ids, drop = FALSE]
  # Clear row/colnames to avoid BGLR confusion with duplicated names
  rownames(K_obs) <- NULL
  colnames(K_obs) <- NULL
  return(K_obs)
}

# ── Build fixed-effect design matrix (FIXED) ─────────────────────────────────
#
# Mirrors the ASReml parameterisation:
#   fixed = trait ~ location
#   random = ~ ... + location:year
#
# Since BGLR has no random effects, we fit location as FIXED and
# year nested within location as FIXED.  BGLR always fits its own
# intercept (mu), so we drop the intercept column from model.matrix().
#
# Result: location main effect  (n_loc - 1 columns)
#       + year:location interaction (n_years_per_loc - 1 columns per location)
#
# This is NOT identical to the ASReml model (where location:year is random
# and therefore shrunk toward zero), but it is the closest BGLR equivalent
# and much more comparable than a single location_year factor.
# Identical to BayesC_utils.R.

build_fixed_design <- function(location_vec, year_vec) {
  loc_factor  <- factor(location_vec)
  year_factor <- factor(year_vec)

  # location main effect (reference = first level, absorbed into BGLR mu)
  Z_loc <- model.matrix(~ loc_factor)[, -1, drop = FALSE]
  colnames(Z_loc) <- paste0("loc_", levels(loc_factor)[-1])

  # year nested within location
  Z_year_in_loc <- NULL
  for (loc in levels(loc_factor)) {
    rows_this_loc <- which(location_vec == loc)
    years_this_loc <- sort(unique(year_vec[rows_this_loc]))
    if (length(years_this_loc) <= 1) next  # nothing to estimate

    # Dummy for each non-reference year within this location
    ref_year <- years_this_loc[1]
    for (yr in years_this_loc[-1]) {
      col_name <- paste0(loc, "_yr_", yr)
      new_col <- integer(length(location_vec))
      new_col[location_vec == loc & year_vec == yr] <- 1L
      if (is.null(Z_year_in_loc)) {
        Z_year_in_loc <- matrix(new_col, ncol = 1,
                                dimnames = list(NULL, col_name))
      } else {
        Z_year_in_loc <- cbind(Z_year_in_loc,
                               setNames(data.frame(new_col), col_name))
      }
    }
  }

  if (!is.null(Z_year_in_loc)) {
    Z_year_in_loc <- as.matrix(Z_year_in_loc)
    Z <- cbind(Z_loc, Z_year_in_loc)
  } else {
    Z <- Z_loc
  }

  return(Z)
}

#' Build environment group vector for heterogeneous residual variances.
#' Identical to BayesC_utils.R.

build_groups <- function(location_year_vec) {
  as.integer(factor(location_year_vec))
}


# ============================================================================
# BGLR RKHS MODEL FITTING
# ============================================================================

#' Fit multi-kernel RKHS model and return observation-level predictions.
#'
#' ETA structure:
#'   1. Location-year fixed effects (contrast-coded design matrix)
#'   2. RKHS kernel at h1 (narrow bandwidth)
#'   3. RKHS kernel at h2 (medium bandwidth)
#'   4. RKHS kernel at h3 (wide bandwidth)
#'
#' NOTE: Unlike BayesC, BGLR does not support heterogeneous residual
#' variances (groups) with RKHS model type. A single residual variance
#' is estimated across all environments.
#'
#' @param model_data Data frame with sample.id, location_year, location,
#'        and the trait column (test cells already set to NA).
#' @param trait Character: trait column name.
#' @param K_geno_list List of 3 genotype-level kernel matrices.
#' @param nIter,burnIn,thin MCMC settings.
#' @param saveAt Character: prefix for BGLR temp files.
#' @return Data frame with sample.id, location_year, location, predicted;
#'         or NULL on failure.

fit_rkhs_and_predict <- function(model_data, trait, K_geno_list,
                                  nIter, burnIn, thin, saveAt) {

  non_na_count <- sum(!is.na(model_data[[trait]]))
  if (non_na_count < 50) {
    warning(paste("Insufficient training data for", trait, ":", non_na_count, "obs"))
    return(NULL)
  }

  # Build observation-level kernels
  sample_ids <- model_data$sample.id
  K_obs_list <- lapply(K_geno_list, function(K) {
    build_obs_kernel(sample_ids, K)
  })

  # Build fixed-effect design matrix
  Z_fix  <- build_fixed_design(model_data$location, model_data$year)

  # Response vector
  y <- model_data[[trait]]

  # ETA: fixed effects (location + year-within-location) + 3 RKHS kernels
  ETA <- list(
    list(X = Z_fix,         model = "FIXED"),
    list(K = K_obs_list[[1]], model = "RKHS"),
    list(K = K_obs_list[[2]], model = "RKHS"),
    list(K = K_obs_list[[3]], model = "RKHS")
  )

  # NOTE: BGLR does not support heterogeneous residual variances (groups)
  # with RKHS model type. Unlike BayesC, RKHS uses a single residual
  # variance across all environments.
  fm <- tryCatch({
    BGLR(
      y       = y,
      ETA     = ETA,
      nIter   = nIter,
      burnIn  = burnIn,
      thin    = thin,
      verbose = FALSE,
      saveAt  = saveAt
    )
  }, error = function(e) {
    warning(paste("BGLR fitting error:", e$message))
    return(NULL)
  })

  if (is.null(fm)) return(NULL)

  # Optionally report kernel variance components
  varU <- sapply(fm$ETA[2:4], function(x) x$varU)
  cat("      varU (h1/h2/h3):", paste(round(varU, 3), collapse = " / "), "\n")

  pred_df <- data.frame(
    sample.id     = model_data$sample.id,
    location_year = model_data$location_year,
    location      = model_data$location,
    predicted     = fm$yHat,
    stringsAsFactors = FALSE
  )

  # Clean up BGLR temp files
  bglr_files <- list.files(
    pattern = paste0("^", gsub("([.|()\\^{}+$*?])", "\\\\\\1", saveAt)),
    full.names = TRUE)
  if (length(bglr_files) > 0) file.remove(bglr_files)

  return(pred_df)
}


# ============================================================================
# EVALUATION METRICS (identical to BayesC_utils.R)
# ============================================================================

calculate_ndcg <- function(y_true, y_pred, k = 10, lower_is_better = FALSE) {
  k <- min(k, length(y_true))

  if (lower_is_better) {
    y_true <- -y_true
    y_pred <- -y_pred
  }

  min_val <- min(c(y_true, y_pred))
  if (min_val < 0) {
    y_true_pos <- y_true - min_val + 1e-6
    y_pred_pos <- y_pred - min_val + 1e-6
  } else {
    y_true_pos <- y_true
    y_pred_pos <- y_pred
  }

  pred_order <- order(y_pred_pos, decreasing = TRUE)
  relevance <- y_true_pos[pred_order[1:k]]
  dcg <- sum(relevance / log2(seq_len(k) + 1))

  ideal_relevance <- sort(y_true_pos, decreasing = TRUE)[1:k]
  idcg <- sum(ideal_relevance / log2(seq_len(k) + 1))

  if (idcg == 0) return(0)
  return(dcg / idcg)
}

evaluate_predictions <- function(y_true, y_pred, trait_name = NULL) {
  if (length(y_true) < 2) {
    return(list(pearson = NA, spearman = NA, ndcg_at_10 = NA))
  }

  lower_is_better <- if (!is.null(trait_name)) {
    trait_name %in% LOWER_IS_BETTER_TRAITS
  } else FALSE

  tryCatch({
    list(
      pearson    = cor(y_true, y_pred, use = "complete.obs"),
      spearman   = cor(y_true, y_pred, use = "complete.obs", method = "spearman"),
      ndcg_at_10 = calculate_ndcg(y_true, y_pred, k = 10,
                                   lower_is_better = lower_is_better)
    )
  }, error = function(e) {
    list(pearson = NA, spearman = NA, ndcg_at_10 = NA)
  })
}


# ============================================================================
# PER-LOCATION-YEAR EVALUATION (identical to BayesC_utils.R)
# ============================================================================

evaluate_per_location_year <- function(observed_data, pred_values, test_indices,
                                        trait, min_genotypes = 10) {

  results <- data.frame(
    trait = character(), location_year = character(), location = character(),
    pearson = numeric(), spearman = numeric(), ndcg_at_10 = numeric(),
    n_test_genotypes = integer(), stringsAsFactors = FALSE
  )

  obs_df <- data.frame(
    sample.id     = observed_data$sample.id[test_indices],
    location_year = observed_data$location_year[test_indices],
    location      = observed_data$location[test_indices],
    observed      = observed_data[[trait]][test_indices],
    stringsAsFactors = FALSE
  )

  pred_df <- pred_values[test_indices, , drop = FALSE]

  merged <- data.frame(
    sample.id     = obs_df$sample.id,
    location_year = obs_df$location_year,
    location      = obs_df$location,
    observed      = obs_df$observed,
    predicted     = pred_df$predicted,
    stringsAsFactors = FALSE
  )
  merged <- merged[!is.na(merged$observed) & !is.na(merged$predicted), ]

  if (nrow(merged) == 0) return(list(metrics = results, predictions = merged))

  for (ly in unique(merged$location_year)) {
    ly_data <- merged[merged$location_year == ly, ]
    n_geno <- length(unique(ly_data$sample.id))

    if (n_geno < min_genotypes) {
      cat("      Skipping", ly, "- only", n_geno, "test genotypes",
          "(minimum:", min_genotypes, ")\n")
      next
    }

    eval_res <- evaluate_predictions(ly_data$observed, ly_data$predicted,
                                      trait_name = trait)

    if (!is.na(eval_res$pearson)) {
      results <- rbind(results, data.frame(
        trait = trait, location_year = ly, location = ly_data$location[1],
        pearson = eval_res$pearson, spearman = eval_res$spearman,
        ndcg_at_10 = eval_res$ndcg_at_10,
        n_test_genotypes = n_geno, stringsAsFactors = FALSE
      ))
      cat("      ", ly, "- r:", round(eval_res$pearson, 3),
          "| rho:", round(eval_res$spearman, 3),
          "| NDCG@10:", round(eval_res$ndcg_at_10, 3),
          "| N:", n_geno, "\n")
    }
  }

  return(list(metrics = results, predictions = merged))
}

cat("RKHS_utils.R loaded successfully.\n")