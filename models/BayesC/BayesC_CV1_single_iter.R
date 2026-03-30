################################################################################
# BayesC_CV1_single_iter.R
#
# Runs ONE iteration of CV1 (new genotypes in known location-years) for
# one trait. Loops over all k folds within this iteration.
#
# Usage:
#   Rscript BayesC_CV1_single_iter.R <trait> <iteration> [marker_file]
#   e.g.  Rscript BayesC_CV1_single_iter.R DTF_blue 3
#   e.g.  Rscript BayesC_CV1_single_iter.R DTF_blue 3 pruned05_AUSPAK_for_bayesC.raw
#
# Output:
#   cv_results_CV1_<trait>_iter<NN>_BayesC.csv
#
# The iteration number determines the random seed (1000 + iteration),
# so results are fully reproducible and identical to the monolithic script.
################################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript BayesC_CV1_single_iter.R <trait> <iteration> [marker_file]")
}

TRAIT       <- args[1]
ITERATION   <- as.integer(args[2])
MARKER_FILE <- if (length(args) >= 3) args[3] else "pruned05_AUSPAK_for_bayesC.raw"

if (is.na(ITERATION) || ITERATION < 1) {
  stop("Iteration must be a positive integer, got: '", args[2], "'")
}

cat("================================================================\n")
cat("BayesC CV1 — Single Iteration\n")
cat("Trait:", TRAIT, "| Iteration:", ITERATION, "\n")
cat("Marker file:", MARKER_FILE, "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

# Source shared utilities (must be in the same directory or working dir)
source("BayesC_utils.R")

# ── Load and prepare data ────────────────────────────────────────────────────
dat <- load_and_prepare_data(TRAIT, MARKER_FILE)
pheno_scaled <- dat$pheno
X_geno       <- dat$X_geno

# ── CV parameters ────────────────────────────────────────────────────────────
k_folds       <- 5
min_genotypes <- 10

# MCMC settings — calibrate from diagnostic run
# (override via environment variables for testing)
nIter  <- as.integer(Sys.getenv("BAYESC_NITER",  unset = "15000"))
burnIn <- as.integer(Sys.getenv("BAYESC_BURNIN", unset = "5000"))
thin   <- as.integer(Sys.getenv("BAYESC_THIN",   unset = "5"))

# ── Run single iteration of CV1 ─────────────────────────────────────────────

pheno_observed <- pheno_scaled  # unmasked copy
genotypes      <- unique(pheno_scaled$sample.id)
n_genotypes    <- length(genotypes)

cat("Data:", n_genotypes, "genotypes,",
    length(unique(pheno_scaled$location_year)), "location-years,",
    nrow(pheno_scaled), "observations\n")
cat("Running iteration", ITERATION, "with", k_folds, "folds\n\n")

current_seed <- 1000 + ITERATION
set.seed(current_seed)
cat("Seed:", current_seed, "\n")

fold_assignments <- sample(rep(1:k_folds, length.out = n_genotypes))
names(fold_assignments) <- genotypes

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
  
  test_genotypes <- names(fold_assignments)[fold_assignments == fold]
  
  tryCatch({
    model_data <- pheno_scaled
    mask_rows  <- which(model_data$sample.id %in% test_genotypes)
    model_data[[TRAIT]][mask_rows] <- NA
    
    if (sum(!is.na(model_data[[TRAIT]])) < 50) {
      warning("Too few training observations"); next
    }
    
    saveAt <- paste0("BayesC_CV1_", TRAIT, "_i", ITERATION, "_f", fold, "_")
    
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
      ly_results$cv_scheme <- "CV1"
      cv_results <- rbind(cv_results, ly_results)
    }

    if (nrow(eval_out$predictions) > 0) {
      preds <- eval_out$predictions
      preds$trait     <- TRAIT
      preds$iteration <- ITERATION
      preds$fold      <- fold
      preds$seed      <- current_seed
      preds$cv_scheme <- "CV1"
      cv_predictions <- rbind(cv_predictions, preds)
    }
    
  }, error = function(e) {
    warning(paste("Error CV1 iter", ITERATION, "fold", fold, ":", e$message))
  })
}

# ── Save results ─────────────────────────────────────────────────────────────

# Zero-pad iteration number for clean file sorting
outfile <- sprintf("cv_results_CV1_%s_iter%02d_BayesC.csv", TRAIT, ITERATION)
write.csv(cv_results, outfile, row.names = FALSE)

pred_outfile <- sprintf("predictions_CV1_%s_iter%02d_BayesC.csv", TRAIT, ITERATION)
write.csv(cv_predictions, pred_outfile, row.names = FALSE)

cat("\n================================================================\n")
cat("CV1 iteration", ITERATION, "complete for trait:", TRAIT, "\n")
cat("Evaluations:", nrow(cv_results), "\n")
cat("Results saved to:", outfile, "\n")
cat("Predictions saved to:", pred_outfile, "(", nrow(cv_predictions), "rows )\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n")