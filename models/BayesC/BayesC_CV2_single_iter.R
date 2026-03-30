################################################################################
# BayesC_CV2_single_iter.R
#
# Runs ONE iteration of CV2 (sparse testing — random observation-level cells
# masked) for one trait. Uses 5-fold CV at the observation level: observed
# cells are split into 5 folds, one fold is masked per fit, and predictions
# are evaluated per location-year within each fold.
#
# This produces the same number of evaluation rows per iteration as CV1
# (5 folds × location-years), enabling direct comparison of CV schemes.
#
# Usage:
#   Rscript BayesC_CV2_single_iter.R <trait> <iteration> [marker_file]
#   e.g.  Rscript BayesC_CV2_single_iter.R DTF_blue 3
#
# Output:
#   cv_results_CV2_<trait>_iter<NN>_BayesC.csv
#
# The iteration number determines the random seed (2000 + iteration),
# identical to the monolithic script's seed scheme.
################################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript BayesC_CV2_single_iter.R <trait> <iteration> [marker_file]")
}

TRAIT       <- args[1]
ITERATION   <- as.integer(args[2])
MARKER_FILE <- if (length(args) >= 3) args[3] else "pruned05_AUSPAK_for_bayesC.raw"

if (is.na(ITERATION) || ITERATION < 1) {
  stop("Iteration must be a positive integer, got: '", args[2], "'")
}

cat("================================================================\n")
cat("BayesC CV2 — Single Iteration (Sparse Testing, 5-fold)\n")
cat("Trait:", TRAIT, "| Iteration:", ITERATION, "\n")
cat("Marker file:", MARKER_FILE, "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

source("BayesC_utils.R")

# ── Load and prepare data ────────────────────────────────────────────────────
dat <- load_and_prepare_data(TRAIT, MARKER_FILE)
pheno_scaled <- dat$pheno
X_geno       <- dat$X_geno

# ── CV parameters ────────────────────────────────────────────────────────────
k_folds       <- 5
min_genotypes <- 10

# MCMC settings (override via environment variables for testing)
nIter  <- as.integer(Sys.getenv("BAYESC_NITER",  unset = "15000"))
burnIn <- as.integer(Sys.getenv("BAYESC_BURNIN", unset = "5000"))
thin   <- as.integer(Sys.getenv("BAYESC_THIN",   unset = "5"))

# ── Run single iteration of CV2 ─────────────────────────────────────────────

pheno_observed <- pheno_scaled  # unmasked copy

cat("Data:", length(unique(pheno_scaled$sample.id)), "genotypes,",
    length(unique(pheno_scaled$location_year)), "location-years,",
    nrow(pheno_scaled), "observations\n")
cat("Running iteration", ITERATION, "with", k_folds, "folds\n\n")

current_seed <- 2000 + ITERATION
set.seed(current_seed)
cat("Seed:", current_seed, "\n")

# Identify observed (non-NA) indices
observed_idx <- which(!is.na(pheno_scaled[[TRAIT]]))
n_observed   <- length(observed_idx)

if (n_observed < 100) {
  stop(paste("Too few observations for", TRAIT, ":", n_observed))
}

# Stratified fold assignment: assign observed cells to folds within each
# location-year so that each fold has balanced representation across
# environments. This prevents a fold from accidentally containing most
# observations from a small location-year.
fold_assignments <- integer(nrow(pheno_scaled))  # 0 for unobserved

location_years_obs <- pheno_scaled$location_year[observed_idx]

for (ly in unique(location_years_obs)) {
  ly_obs_idx <- observed_idx[location_years_obs == ly]
  n_ly <- length(ly_obs_idx)
  # Assign folds within this location-year
  fold_assignments[ly_obs_idx] <- sample(rep(1:k_folds, length.out = n_ly))
}

# Report fold sizes per location-year
cat("\nFold sizes per location-year:\n")
fold_table <- table(
  pheno_scaled$location_year[observed_idx],
  fold_assignments[observed_idx]
)
print(fold_table)
cat("\n")

cv_results <- data.frame(
  iteration = integer(), fold = integer(), trait = character(),
  location_year = character(), location = character(),
  pearson = numeric(), spearman = numeric(), ndcg_at_10 = numeric(),
  seed = integer(), n_test_genotypes = integer(), cv_scheme = character(),
  stringsAsFactors = FALSE
)

cv_predictions <- data.frame(
  sample.id = character(), location_year = character(),
  location = character(), observed = numeric(), predicted = numeric(),
  trait = character(), iteration = integer(), fold = integer(),
  seed = integer(), cv_scheme = character(),
  stringsAsFactors = FALSE
)

for (fold in 1:k_folds) {
  cat("\n  Fold", fold, "of", k_folds, "\n")
  
  mask_rows <- which(fold_assignments == fold)
  
  tryCatch({
    model_data <- pheno_scaled
    model_data[[TRAIT]][mask_rows] <- NA
    
    n_training <- sum(!is.na(model_data[[TRAIT]]))
    if (n_training < 50) {
      warning("Too few training observations: ", n_training); next
    }
    
    # Safety check: ensure each location-year retains enough training data
    ly_train_counts <- table(
      pheno_scaled$location_year[which(!is.na(model_data[[TRAIT]]))]
    )
    thin_lys <- names(ly_train_counts)[ly_train_counts < 10]
    if (length(thin_lys) > 0) {
      cat("    Warning: location-years with <10 training obs:",
          paste(thin_lys, collapse = ", "), "\n")
    }
    
    saveAt <- paste0("BayesC_CV2_", TRAIT, "_i", ITERATION, "_f", fold, "_")
    
    pred_values <- fit_bayesc_and_predict(
      model_data = model_data, trait = TRAIT, X_geno = X_geno,
      nIter = nIter, burnIn = burnIn, thin = thin, saveAt = saveAt
    )
    
    if (is.null(pred_values)) next
    
    test_indices <- mask_rows
    
    eval_out <- evaluate_per_location_year(
      observed_data = pheno_observed,
      pred_values   = pred_values,
      test_indices  = test_indices,
      trait = TRAIT, min_genotypes = min_genotypes
    )

    ly_results <- eval_out$metrics
    if (nrow(ly_results) > 0) {
      ly_results$iteration <- ITERATION
      ly_results$fold      <- fold
      ly_results$seed      <- current_seed
      ly_results$cv_scheme <- "CV2"
      cv_results <- rbind(cv_results, ly_results)
    }

    if (nrow(eval_out$predictions) > 0) {
      preds <- eval_out$predictions
      preds$trait     <- TRAIT
      preds$iteration <- ITERATION
      preds$fold      <- fold
      preds$seed      <- current_seed
      preds$cv_scheme <- "CV2"
      cv_predictions <- rbind(cv_predictions, preds)
    }
    
  }, error = function(e) {
    warning(paste("Error CV2 iter", ITERATION, "fold", fold, ":", e$message))
  })
}

# ── Save results ─────────────────────────────────────────────────────────────

outfile <- sprintf("cv_results_CV2_%s_iter%02d_BayesC.csv", TRAIT, ITERATION)
write.csv(cv_results, outfile, row.names = FALSE)

pred_outfile <- sprintf("predictions_CV2_%s_iter%02d_BayesC.csv", TRAIT, ITERATION)
write.csv(cv_predictions, pred_outfile, row.names = FALSE)

cat("\n================================================================\n")
cat("CV2 iteration", ITERATION, "complete for trait:", TRAIT, "\n")
cat("Evaluations:", nrow(cv_results), "\n")
cat("Results saved to:", outfile, "\n")
cat("Predictions saved to:", pred_outfile, "(", nrow(cv_predictions), "rows )\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n")