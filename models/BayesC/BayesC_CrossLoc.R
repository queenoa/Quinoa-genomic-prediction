################################################################################
# BayesC_CrossLoc.R
#
# Runs cross-location transferability for one trait.
# Trains on each location separately, predicts into all other location-years.
# Results and predictions are saved incrementally after each training location.
#
# Usage:
#   Rscript BayesC_CrossLoc.R <trait> [marker_file]
#   e.g.  Rscript BayesC_CrossLoc.R DTF_blue
#
# Output:
#   cv_results_CrossLoc_<trait>_BayesC.csv
#   predictions_CrossLoc_<trait>_BayesC.csv
################################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript BayesC_CrossLoc.R <trait> [marker_file]")
}

TRAIT       <- args[1]
MARKER_FILE <- if (length(args) >= 2) args[2] else "pruned05_AUSPAK_for_bayesC.raw"

cat("================================================================\n")
cat("BayesC Cross-Location Transferability\n")
cat("Trait:", TRAIT, "\n")
cat("Marker file:", MARKER_FILE, "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

source("BayesC_utils.R")

# ── Load and prepare data ────────────────────────────────────────────────────
dat <- load_and_prepare_data(TRAIT, MARKER_FILE)
pheno_scaled <- dat$pheno
X_geno       <- dat$X_geno

# ── Parameters ───────────────────────────────────────────────────────────────
min_genotypes <- 10
nIter  <- as.integer(Sys.getenv("BAYESC_NITER",  unset = "15000"))
burnIn <- as.integer(Sys.getenv("BAYESC_BURNIN", unset = "5000"))
thin   <- as.integer(Sys.getenv("BAYESC_THIN",   unset = "5"))

# ── Output file paths ───────────────────────────────────────────────────────
outfile_results <- paste0("cv_results_CrossLoc_", TRAIT, "_BayesC.csv")
outfile_preds   <- paste0("predictions_CrossLoc_", TRAIT, "_BayesC.csv")

# ── Cross-location transferability ──────────────────────────────────────────

locations <- sort(unique(pheno_scaled$location))
cat("Locations:", paste(locations, collapse = ", "), "\n\n")

# Write CSV headers (overwrite any previous partial run)
results_header <- data.frame(
  iteration = integer(), fold = integer(), trait = character(),
  location_year = character(), location = character(),
  pearson = numeric(), spearman = numeric(), ndcg_at_10 = numeric(),
  seed = integer(), n_test_genotypes = integer(), cv_scheme = character(),
  stringsAsFactors = FALSE
)
write.csv(results_header, outfile_results, row.names = FALSE)

preds_header <- data.frame(
  sample.id = character(), location_year = character(),
  location = character(), observed = numeric(), predicted = numeric(),
  trait = character(), train_location = character(), cv_scheme = character(),
  stringsAsFactors = FALSE
)
write.csv(preds_header, outfile_preds, row.names = FALSE)

n_evaluations <- 0

for (train_loc in locations) {
  predict_locs <- setdiff(locations, train_loc)
  cat("Training on:", train_loc, "-> Predicting:",
      paste(predict_locs, collapse = ", "), "\n")

  train_data <- pheno_scaled[pheno_scaled$location == train_loc, ]

  n_train_obs <- sum(!is.na(train_data[[TRAIT]]))
  if (n_train_obs < 50) {
    cat("  Skipping - only", n_train_obs, "training observations\n")
    next
  }

  tryCatch({
    X_train <- build_obs_marker_matrix(train_data$sample.id, X_geno)

    train_years <- unique(train_data$year)

    if (length(train_years) > 1) {
      Z_year <- model.matrix(~ factor(train_data$year))[, -1, drop = FALSE]
      ETA <- list(
        list(X = Z_year,  model = "FIXED"),
        list(X = X_train, model = "BayesC")
      )
    } else {
      ETA <- list(
        list(X = X_train, model = "BayesC")
      )
    }

    y_train <- train_data[[TRAIT]]

    if (length(train_years) > 1) {
      groups_train <- as.integer(factor(train_data$location_year))
    } else {
      groups_train <- NULL
    }

    saveAt <- paste0("BayesC_CrossLoc_", TRAIT, "_train", train_loc, "_")

    fm <- BGLR(
      y       = y_train,
      ETA     = ETA,
      groups  = groups_train,
      nIter   = nIter,
      burnIn  = burnIn,
      thin    = thin,
      verbose = FALSE,
      saveAt  = saveAt
    )

    bayesc_idx <- length(ETA)
    beta_hat <- fm$ETA[[bayesc_idx]]$b
    cat("  Extracted marker effects (", length(beta_hat), " SNPs)\n")

    # Clean up BGLR temp files
    bglr_files <- list.files(
      pattern = paste0("^", gsub("([.|()\\^{}+$*?])", "\\\\\\1", saveAt)),
      full.names = TRUE)
    if (length(bglr_files) > 0) file.remove(bglr_files)

    # Predict into target location-years
    target_lys <- unique(pheno_scaled$location_year[
      pheno_scaled$location %in% predict_locs])

    preds_this_loc <- data.frame()
    results_this_loc <- data.frame()

    for (ly in target_lys) {
      target_data <- pheno_scaled[pheno_scaled$location_year == ly, ]
      target_obs  <- target_data[!is.na(target_data[[TRAIT]]), ]

      if (nrow(target_obs) == 0) next

      X_target <- X_geno[target_obs$sample.id, , drop = FALSE]
      gebv <- as.numeric(fm$mu + X_target %*% beta_hat)

      cv_label <- paste0("CrossLoc_", train_loc, "->",
                         paste(predict_locs, collapse = "+"))

      # Collect predictions
      preds_this_loc <- rbind(preds_this_loc, data.frame(
        sample.id      = target_obs$sample.id,
        location_year  = target_obs$location_year,
        location       = target_obs$location,
        observed       = target_obs[[TRAIT]],
        predicted      = gebv,
        trait           = TRAIT,
        train_location  = train_loc,
        cv_scheme       = cv_label,
        stringsAsFactors = FALSE
      ))

      n_geno <- length(unique(target_obs$sample.id))

      if (n_geno < min_genotypes) {
        cat("    Skipping", ly, "- only", n_geno, "genotypes\n")
        next
      }

      eval_res <- evaluate_predictions(
        target_obs[[TRAIT]], gebv, trait_name = TRAIT
      )

      if (!is.na(eval_res$pearson)) {
        target_loc <- unique(target_obs$location)
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
      write.table(preds_this_loc, outfile_preds, append = TRUE, sep = ",",
                  row.names = FALSE, col.names = FALSE, quote = TRUE)
      cat("  Predictions appended (", nrow(preds_this_loc), "rows )\n")
    }
    if (nrow(results_this_loc) > 0) {
      write.table(results_this_loc, outfile_results, append = TRUE, sep = ",",
                  row.names = FALSE, col.names = FALSE, quote = TRUE)
      n_evaluations <- n_evaluations + nrow(results_this_loc)
    }

  }, error = function(e) {
    warning(paste("Error cross-location for train_loc", train_loc, ":", e$message))
  })
}

cat("\n================================================================\n")
cat("CrossLoc complete for trait:", TRAIT, "\n")
cat("CrossLoc evaluations:", n_evaluations, "\n")
cat("Results:", outfile_results, "\n")
cat("Predictions:", outfile_preds, "\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n")
