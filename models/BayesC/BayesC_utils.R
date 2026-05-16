################################################################################
# BayesC_utils.R — Shared utility functions for BayesC multi-environment CV
#
# Sourced by:
#   BayesC_CV1_single_iter.R
#   BayesC_CV2_single_iter.R
#   BayesC_CV0_CrossLoc.R
#
# Contains:
#   - Data loading and preprocessing
#   - Z-score standardisation
#   - Marker matrix construction
#   - BGLR model fitting wrapper
#   - Evaluation metrics (Pearson, Spearman, NDCG@10)
#   - Summarisation helpers
################################################################################

library(BGLR)
library(dplyr)
library(data.table)

# ── Global constants ─────────────────────────────────────────────────────────

VALID_TRAITS <- c('DTF', 'DTH', 'PtHt', 'PcleLng',
                  'SdLen', 'TGW', 'SdW_z')

LOWER_IS_BETTER_TRAITS <- c('DTF', 'DTH', 'PtHt')

# ── Data loading and preprocessing ───────────────────────────────────────────

load_and_prepare_data <- function(trait, marker_file, pheno_file = "../AUSPAK_phenotypes_GP_input.csv") {
  
  # Validate trait
  if (!trait %in% VALID_TRAITS) {
    stop("Invalid trait: '", trait, "'\n",
         "  Valid traits: ", paste(VALID_TRAITS, collapse = ", "))
  }
  
  if (!file.exists(marker_file)) {
    stop("Marker file not found: ", marker_file)
  }
  if (!file.exists(pheno_file)) {
    stop("Phenotype file not found: ", pheno_file)
  }
  
  # Load phenotype data
  pheno <- read.csv(pheno_file, stringsAsFactors = FALSE)
  
  required_cols <- c("sample.id", "location", "year", "location_year", trait)
  missing_cols <- setdiff(required_cols, names(pheno))
  if (length(missing_cols) > 0) {
    stop("Missing columns in phenotype file: ", paste(missing_cols, collapse = ", "))
  }
  
  cat("Phenotype data:", nrow(pheno), "rows,", ncol(pheno), "columns\n")
  cat("Location-years:", paste(sort(unique(pheno$location_year)), collapse = ", "), "\n")
  
  # Load marker matrix
  raw <- fread(marker_file, header = TRUE, check.names = FALSE)
  geno_ids <- as.character(raw$IID)
  X_geno <- as.matrix(raw[, -(1:6)])
  rownames(X_geno) <- geno_ids
  
  cat("Marker matrix:", nrow(X_geno), "genotypes x", ncol(X_geno), "SNPs\n")
  
  # Mean-impute missing genotypes
  na_count <- sum(is.na(X_geno))
  cat("Missing genotypes:", na_count, "(",
      round(mean(is.na(X_geno)) * 100, 4), "% of all marker calls)\n")
  
  if (na_count > 0) {
    imp_vals <- colMeans(X_geno, na.rm = TRUE)
    ix <- which(is.na(X_geno), arr.ind = TRUE)
    X_geno[ix] <- imp_vals[ix[, 2]]
    cat("Mean-imputed", nrow(ix), "missing calls across",
        length(unique(ix[, 2])), "markers.\n")
  }
  stopifnot(sum(is.na(X_geno)) == 0)
  
  rm(raw); gc()

  # centre the marker matrix

  X_geno <- scale(X_geno, center = TRUE, scale = FALSE)
  
  # Align phenotype and genotype data
  pheno$sample.id     <- as.character(pheno$sample.id)
  pheno$location      <- as.character(pheno$location)
  pheno$year          <- as.character(pheno$year)
  pheno$location_year <- as.character(pheno$location_year)
  
  pheno <- pheno[pheno$sample.id %in% geno_ids, ]
  cat("Observations with genotype data:", nrow(pheno), "\n")
  cat("Unique genotypes:", length(unique(pheno$sample.id)), "\n")
  
  stopifnot(all(pheno$sample.id %in% rownames(X_geno)))
  
  # Z-score standardisation by location-year
  cat("Applying z-score standardisation by location-year...\n")
  pheno <- apply_location_year_scaling(pheno, trait)
  
  # Report missingness
  n_obs <- sum(!is.na(pheno[[trait]]))
  n_total <- nrow(pheno)
  cat("Trait", trait, ":", n_obs, "observed values out of", n_total,
      "(", round(100 * (n_total - n_obs) / n_total, 1), "% missing)\n")
  
  # Drop location-years with zero observed values
  ly_summary <- pheno %>%
    group_by(location_year) %>%
    summarise(n_total = n(),
              n_observed = sum(!is.na(.data[[trait]])),
              .groups = 'drop')
  
  empty_lys <- ly_summary$location_year[ly_summary$n_observed == 0]
  if (length(empty_lys) > 0) {
    cat("Dropping location-years with 0 observed", trait, "values:\n")
    cat("  ", paste(empty_lys, collapse = ", "), "\n")
    n_before <- nrow(pheno)
    pheno <- pheno[!pheno$location_year %in% empty_lys, ]
    cat("  Rows:", n_before, "->", nrow(pheno),
        "(removed", n_before - nrow(pheno), "rows)\n")
  }
  
  cat("Remaining location-years:", paste(sort(unique(pheno$location_year)),
                                          collapse = ", "), "\n\n")
  
  return(list(pheno = pheno, X_geno = X_geno))
}

# ── Z-score standardisation by location-year ─────────────────────────────────

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

# ── Build observation-level marker matrix ────────────────────────────────────

build_obs_marker_matrix <- function(sample_ids, X_geno) {
  X_obs <- X_geno[sample_ids, , drop = FALSE]
  rownames(X_obs) <- NULL
  return(X_obs)
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

build_fixed_design <- function(location_vec, year_vec) {
  loc_factor  <- factor(location_vec)
  year_factor <- factor(year_vec)

  # location main effect (reference = first level, absorbed into BGLR mu)
  Z_loc <- model.matrix(~ loc_factor)[, -1, drop = FALSE]
  colnames(Z_loc) <- paste0("loc_", levels(loc_factor)[-1])

  # year nested within location (interaction dummies, drop intercept)
  Z_loc_year <- model.matrix(~ loc_factor:year_factor)[, -1, drop = FALSE]

  # Remove columns that duplicate the location main effect
  # (model.matrix for interaction includes main-effect-like columns)
  # Keep only columns with real interaction variation
  # Safest approach: use the full location:year interaction minus the
  # location main effect, which gives year-within-location contrasts.
  #
  # Equivalent to: model.matrix(~ location/year - 1) then dropping one
  # reference year per location.
  #
  # Simpler and more robust approach: build it manually per location.
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

# ── Backward-compatible wrapper (used nowhere now, kept for safety) ──────────

# build_locyear_design <- function(location_year_vec) {
#   warning("build_locyear_design() is deprecated. Use build_fixed_design() instead.")
#   ly_factor <- factor(location_year_vec)
#   Z <- model.matrix(~ ly_factor)[, -1, drop = FALSE]
#   colnames(Z) <- levels(ly_factor)[-1]
#   return(Z)
# }

# ── Build environment group vector ───────────────────────────────────────────

build_groups <- function(location_year_vec) {
  as.integer(factor(location_year_vec))
}

# ── NDCG@k calculation ──────────────────────────────────────────────────────

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

# ── Evaluate predictions ─────────────────────────────────────────────────────

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

# ── Fit BayesC and return predictions ────────────────────────────────────────

fit_bayesc_and_predict <- function(model_data, trait, X_geno,
                                    nIter, burnIn, thin, saveAt) {
  
  non_na_count <- sum(!is.na(model_data[[trait]]))
  if (non_na_count < 50) {
    warning(paste("Insufficient training data for", trait, ":", non_na_count, "obs"))
    return(NULL)
  }
  
  X_obs  <- build_obs_marker_matrix(model_data$sample.id, X_geno)
  Z_fix  <- build_fixed_design(model_data$location, model_data$year)
  groups <- build_groups(model_data$location_year)
  y      <- model_data[[trait]]
  
  ETA <- list(
    list(X = Z_fix, model = "FIXED"),
    list(X = X_obs, model = "BayesC")
  )
  
  fm <- tryCatch({
    BGLR(
      y       = y,
      ETA     = ETA,
      groups  = groups,
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

# ── Evaluate per location-year ───────────────────────────────────────────────

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
  
  if (nrow(merged) == 0) return(results)
  
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

cat("BayesC_utils.R loaded successfully.\n")