################################################################################
# RKHS_CV0_CrossLoc.R
#
# Runs CV0 (leave-one-location-year-out) and Cross-location transferability
# for one trait using multi-kernel RKHS. Both are deterministic single-pass
# schemes.
#
# Cross-location approach:
#   Fit within-location RKHS (intercept + 3 kernels). Target-location
#   genotypes are included in the kernel matrix with y = NA, so BGLR's
#   internal machinery predicts them from kernel similarity to training
#   genotypes. This is the kernel analogue of the BayesC approach
#   (which extracts beta_hat and computes X*beta for target genotypes),
#   but here prediction flows through the kernel rather than explicit
#   marker effects.
#
# Usage:
#   Rscript RKHS_CV0_CrossLoc.R <trait> [marker_file]
#   e.g.  Rscript RKHS_CV0_CrossLoc.R DTF_blue
#
# Output:
#   cv_results_CV0_<trait>_RKHS.csv
#   predictions_CV0_<trait>_RKHS.csv
#   cv_results_CrossLoc_<trait>_RKHS.csv
#   predictions_CrossLoc_<trait>_RKHS.csv
#
# Reference:
#   Pérez & de los Campos (2014) Genetics 198:483-495
################################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript RKHS_CV0_CrossLoc.R <trait> [marker_file]")
}

TRAIT       <- args[1]
MARKER_FILE <- if (length(args) >= 2) args[2] else "auspak_for_rkhs.raw"

cat("================================================================\n")
cat("RKHS CV0 + Cross-Location — Deterministic Schemes (Multi-Kernel)\n")
cat("Trait:", TRAIT, "\n")
cat("Marker file:", MARKER_FILE, "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

source("RKHS_utils.R")

# ── Load and prepare data ────────────────────────────────────────────────────
dat <- load_and_prepare_data(TRAIT, MARKER_FILE)
pheno_scaled <- dat$pheno
K_geno_list  <- dat$K_geno_list

# ── Parameters ───────────────────────────────────────────────────────────────
min_genotypes <- 10
nIter  <- as.integer(Sys.getenv("RKHS_NITER",  unset = "15000"))
burnIn <- as.integer(Sys.getenv("RKHS_BURNIN", unset = "5000"))
thin   <- as.integer(Sys.getenv("RKHS_THIN",   unset = "5"))


# ============================================================================
# CV0: Leave-one-location-year-out
# ============================================================================

cat("\n=============================================\n")
cat("CV0: Leave-one-location-year-out\n")
cat("=============================================\n\n")

pheno_observed <- pheno_scaled
location_years <- sort(unique(pheno_scaled$location_year))

cat("Data:", length(unique(pheno_scaled$sample.id)), "genotypes,",
    length(location_years), "location-years,",
    nrow(pheno_scaled), "observations\n\n")

# Output file paths
outfile_cv0       <- paste0("cv_results_CV0_", TRAIT, "_RKHS.csv")
outfile_cv0_preds <- paste0("predictions_CV0_", TRAIT, "_RKHS.csv")

# Write CSV headers (overwrite any previous partial run)
cv0_results <- data.frame(
  iteration = integer(), fold = integer(), trait = character(),
  location_year = character(), location = character(),
  pearson = numeric(), spearman = numeric(), ndcg_at_10 = numeric(),
  seed = integer(), n_test_genotypes = integer(), cv_scheme = character(),
  stringsAsFactors = FALSE
)
write.csv(cv0_results, outfile_cv0, row.names = FALSE)

cv0_preds_header <- data.frame(
  sample.id = character(), location_year = character(),
  location = character(), observed = numeric(), predicted = numeric(),
  trait = character(), cv_scheme = character(),
  stringsAsFactors = FALSE
)
write.csv(cv0_preds_header, outfile_cv0_preds, row.names = FALSE)

n_cv0_evaluations <- 0

for (held_out_ly in location_years) {
  cat("Holding out:", held_out_ly, "\n")

  held_out_location <- unique(pheno_scaled$location[
    pheno_scaled$location_year == held_out_ly])

  test_genos_in_ly <- unique(pheno_observed$sample.id[
    pheno_observed$location_year == held_out_ly &
      !is.na(pheno_observed[[TRAIT]])])

  if (length(test_genos_in_ly) < min_genotypes) {
    cat("  Skipping - only", length(test_genos_in_ly),
        "genotypes with observed values\n")
    next
  }

  tryCatch({
    model_data <- pheno_scaled
    mask_rows  <- which(model_data$location_year == held_out_ly)
    model_data[[TRAIT]][mask_rows] <- NA

    if (sum(!is.na(model_data[[TRAIT]])) < 50) {
      warning("Too few training observations after masking"); next
    }

    saveAt <- paste0("RKHS_CV0_", TRAIT, "_", held_out_ly, "_")

    pred_values <- fit_rkhs_and_predict(
      model_data  = model_data,
      trait       = TRAIT,
      K_geno_list = K_geno_list,
      nIter = nIter, burnIn = burnIn, thin = thin, saveAt = saveAt
    )

    if (is.null(pred_values)) next

    test_idx <- which(pheno_observed$location_year == held_out_ly &
                        pheno_observed$sample.id %in% test_genos_in_ly)

    obs_vals  <- pheno_observed[[TRAIT]][test_idx]
    pred_vals <- pred_values$predicted[test_idx]

    valid <- !is.na(obs_vals) & !is.na(pred_vals)
    obs_vals  <- obs_vals[valid]
    pred_vals <- pred_vals[valid]
    n_geno    <- length(unique(pheno_observed$sample.id[test_idx][valid]))

    # Append predictions incrementally
    if (sum(valid) > 0) {
      preds_df <- data.frame(
        sample.id     = pheno_observed$sample.id[test_idx][valid],
        location_year = pheno_observed$location_year[test_idx][valid],
        location      = pheno_observed$location[test_idx][valid],
        observed      = obs_vals,
        predicted     = pred_vals,
        trait          = TRAIT,
        cv_scheme      = "CV0",
        stringsAsFactors = FALSE
      )
      write.table(preds_df, outfile_cv0_preds, append = TRUE, sep = ",",
                  row.names = FALSE, col.names = FALSE, quote = TRUE)
      cat("  Predictions appended (", nrow(preds_df), "rows )\n")
    }

    if (n_geno < min_genotypes) {
      cat("  Skipping metrics - only", n_geno, "genotypes with predictions\n")
      next
    }

    eval_res <- evaluate_predictions(obs_vals, pred_vals, trait_name = TRAIT)

    if (!is.na(eval_res$pearson)) {
      results_df <- data.frame(
        iteration = NA, fold = NA, trait = TRAIT,
        location_year = held_out_ly, location = held_out_location,
        pearson = eval_res$pearson, spearman = eval_res$spearman,
        ndcg_at_10 = eval_res$ndcg_at_10,
        seed = NA, n_test_genotypes = n_geno,
        cv_scheme = "CV0", stringsAsFactors = FALSE
      )
      write.table(results_df, outfile_cv0, append = TRUE, sep = ",",
                  row.names = FALSE, col.names = FALSE, quote = TRUE)
      n_cv0_evaluations <- n_cv0_evaluations + 1

      cat("  ", held_out_ly, "- r:", round(eval_res$pearson, 3),
          "| rho:", round(eval_res$spearman, 3),
          "| NDCG@10:", round(eval_res$ndcg_at_10, 3),
          "| N:", n_geno, "\n")
    }

  }, error = function(e) {
    warning(paste("Error CV0 for", held_out_ly, ":", e$message))
  })
}

cat("\nCV0 results saved to:", outfile_cv0, "\n")
cat("CV0 predictions saved to:", outfile_cv0_preds, "\n")
cat("CV0 evaluations:", n_cv0_evaluations, "\n")


# ============================================================================
# Cross-location transferability
# ============================================================================
#
# Strategy: For each training location, build a dataset containing ALL
# genotypes (from both locations) but set y = NA for all rows from the
# target location. Include only the intercept + 3 RKHS kernels (no
# location-year fixed effects, since the target location cannot contribute
# to estimating them — mirrors BayesC cross-location logic).
#
# If the training location has multiple years, a year fixed effect is
# included (contrast-coded) and groups are used for year-specific
# residual variances, matching BayesC.
#
# BGLR predicts for target genotypes via kernel similarity to training
# genotypes (the kernel rows connecting target to training genotypes
# carry the genetic information across locations).

cat("\n=============================================\n")
cat("Cross-location transferability\n")
cat("=============================================\n\n")

locations <- sort(unique(pheno_scaled$location))
cat("Locations:", paste(locations, collapse = ", "), "\n\n")

# Output file paths
outfile_xl       <- paste0("cv_results_CrossLoc_", TRAIT, "_RKHS.csv")
outfile_xl_preds <- paste0("predictions_CrossLoc_", TRAIT, "_RKHS.csv")

# Write CSV headers
crossloc_results <- data.frame(
  iteration = integer(), fold = integer(), trait = character(),
  location_year = character(), location = character(),
  pearson = numeric(), spearman = numeric(), ndcg_at_10 = numeric(),
  seed = integer(), n_test_genotypes = integer(), cv_scheme = character(),
  stringsAsFactors = FALSE
)
write.csv(crossloc_results, outfile_xl, row.names = FALSE)

xl_preds_header <- data.frame(
  sample.id = character(), location_year = character(),
  location = character(), observed = numeric(), predicted = numeric(),
  trait = character(), train_location = character(), cv_scheme = character(),
  stringsAsFactors = FALSE
)
write.csv(xl_preds_header, outfile_xl_preds, row.names = FALSE)

n_xl_evaluations <- 0

for (train_loc in locations) {
  predict_locs <- setdiff(locations, train_loc)
  cat("Training on:", train_loc, "-> Predicting:",
      paste(predict_locs, collapse = ", "), "\n")

  # Build combined dataset: training location (observed) + target location (NA)
  train_data  <- pheno_scaled[pheno_scaled$location == train_loc, ]
  target_data <- pheno_scaled[pheno_scaled$location %in% predict_locs, ]

  n_train_obs <- sum(!is.na(train_data[[TRAIT]]))
  if (n_train_obs < 50) {
    cat("  Skipping - only", n_train_obs, "training observations\n")
    next
  }

  tryCatch({
    # Combine train + target; mask all target phenotypes
    combined_data <- rbind(train_data, target_data)
    target_rows <- which(combined_data$location %in% predict_locs)
    combined_data[[TRAIT]][target_rows] <- NA

    # Build observation-level kernels from the combined dataset
    sample_ids <- combined_data$sample.id
    K_obs_list <- lapply(K_geno_list, function(K) {
      build_obs_kernel(sample_ids, K)
    })

    # Year fixed effect within training location (if >1 year)
    # For RKHS cross-location, we use a simpler model:
    # intercept + kernels only (no location-year fixed effects,
    # since target location-years cannot be estimated).
    # If the training location has multiple years, include year
    # as a fixed effect to account for year-to-year variation
    # within the training data.
    train_years <- unique(train_data$year)

    if (length(train_years) > 1) {
      # Build year design matrix for the COMBINED dataset.
      # Target-location rows get 0s (absorbed into intercept).
      # This is safe because their y is NA — they don't contribute
      # to fixed-effect estimation.
      year_vec <- combined_data$year
      # Map target-location years to the reference level so they
      # get zero columns (BGLR intercept absorbs them)
      ref_year <- levels(factor(train_data$year))[1]
      year_vec[target_rows] <- ref_year
      Z_year <- model.matrix(~ factor(year_vec))[, -1, drop = FALSE]

      ETA <- list(
        list(X = Z_year,          model = "FIXED"),
        list(K = K_obs_list[[1]], model = "RKHS"),
        list(K = K_obs_list[[2]], model = "RKHS"),
        list(K = K_obs_list[[3]], model = "RKHS")
      )
    } else {
      ETA <- list(
        list(K = K_obs_list[[1]], model = "RKHS"),
        list(K = K_obs_list[[2]], model = "RKHS"),
        list(K = K_obs_list[[3]], model = "RKHS")
      )
    }

    y <- combined_data[[TRAIT]]

    # NOTE: BGLR does not support heterogeneous residual variances (groups)
    # with RKHS model type — uses a single residual variance.

    saveAt <- paste0("RKHS_CrossLoc_", TRAIT, "_train", train_loc, "_")

    fm <- BGLR(
      y       = y,
      ETA     = ETA,
      nIter   = nIter,
      burnIn  = burnIn,
      thin    = thin,
      verbose = FALSE,
      saveAt  = saveAt
    )

    # Report kernel variance components
    rkhs_start <- if (length(train_years) > 1) 2 else 1
    varU <- sapply(fm$ETA[rkhs_start:(rkhs_start + 2)], function(x) x$varU)
    cat("  varU (h1/h2/h3):", paste(round(varU, 3), collapse = " / "), "\n")

    # Clean up BGLR temp files
    bglr_files <- list.files(
      pattern = paste0("^", gsub("([.|()\\^{}+$*?])", "\\\\\\1", saveAt)),
      full.names = TRUE)
    if (length(bglr_files) > 0) file.remove(bglr_files)

    # Extract predictions for target-location rows
    pred_all <- fm$yHat

    # Evaluate per target location-year
    target_lys <- unique(combined_data$location_year[target_rows])

    preds_this_loc <- data.frame()
    results_this_loc <- data.frame()

    for (ly in target_lys) {
      ly_rows <- which(combined_data$location_year == ly)

      # Get observed values from the unmasked pheno_scaled
      obs_in_ly <- pheno_scaled[pheno_scaled$location_year == ly &
                                  !is.na(pheno_scaled[[TRAIT]]), ]
      if (nrow(obs_in_ly) == 0) next

      # Match predictions to observed genotypes
      pred_in_ly <- data.frame(
        sample.id = combined_data$sample.id[ly_rows],
        predicted = pred_all[ly_rows],
        stringsAsFactors = FALSE
      )

      merged <- merge(
        obs_in_ly[, c("sample.id", TRAIT, "location")],
        pred_in_ly,
        by = "sample.id"
      )
      merged <- merged[!is.na(merged$predicted), ]

      cv_label <- paste0("CrossLoc_", train_loc, "->",
                         paste(predict_locs, collapse = "+"))

      # Collect predictions
      if (nrow(merged) > 0) {
        preds_this_loc <- rbind(preds_this_loc, data.frame(
          sample.id      = merged$sample.id,
          location_year  = ly,
          location       = merged$location,
          observed       = merged[[TRAIT]],
          predicted      = merged$predicted,
          trait           = TRAIT,
          train_location  = train_loc,
          cv_scheme       = cv_label,
          stringsAsFactors = FALSE
        ))
      }

      n_geno <- length(unique(merged$sample.id))

      if (n_geno < min_genotypes) {
        cat("    Skipping", ly, "- only", n_geno, "genotypes\n")
        next
      }

      eval_res <- evaluate_predictions(
        merged[[TRAIT]], merged$predicted, trait_name = TRAIT
      )

      if (!is.na(eval_res$pearson)) {
        target_loc <- unique(merged$location)
        results_this_loc <- rbind(results_this_loc, data.frame(
          iteration = NA, fold = NA, trait = TRAIT,
          location_year = ly, location = target_loc,
          pearson = eval_res$pearson, spearman = eval_res$spearman,
          ndcg_at_10 = eval_res$ndcg_at_10,
          seed = NA, n_test_genotypes = n_geno,
          cv_scheme = cv_label,
          stringsAsFactors = FALSE
        ))
        cat("    ", ly, "- r:", round(eval_res$pearson, 3),
            "| rho:", round(eval_res$spearman, 3),
            "| NDCG@10:", round(eval_res$ndcg_at_10, 3),
            "| N:", n_geno, "\n")
      }
    }

    # Append incrementally after each training location
    if (nrow(preds_this_loc) > 0) {
      write.table(preds_this_loc, outfile_xl_preds, append = TRUE, sep = ",",
                  row.names = FALSE, col.names = FALSE, quote = TRUE)
      cat("  Predictions appended (", nrow(preds_this_loc), "rows )\n")
    }
    if (nrow(results_this_loc) > 0) {
      write.table(results_this_loc, outfile_xl, append = TRUE, sep = ",",
                  row.names = FALSE, col.names = FALSE, quote = TRUE)
      n_xl_evaluations <- n_xl_evaluations + nrow(results_this_loc)
    }

  }, error = function(e) {
    warning(paste("Error cross-location for train_loc", train_loc, ":", e$message))
  })
}

cat("\nCrossLoc results saved to:", outfile_xl, "\n")
cat("CrossLoc predictions saved to:", outfile_xl_preds, "\n")
cat("CrossLoc evaluations:", n_xl_evaluations, "\n")

# ── Summary ──────────────────────────────────────────────────────────────────

cat("\n================================================================\n")
cat("CV0 + CrossLoc complete for trait:", TRAIT, "\n")
cat("CV0 evaluations:", n_cv0_evaluations, "\n")
cat("CrossLoc evaluations:", n_xl_evaluations, "\n")
cat("Results:", outfile_cv0, ",", outfile_xl, "\n")
cat("Predictions:", outfile_cv0_preds, ",", outfile_xl_preds, "\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n")