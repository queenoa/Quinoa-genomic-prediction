################################################################################
# RKHS_CV1.R
#
# Runs ALL iterations of CV1 (new genotypes in known location-years) for
# one trait using multi-kernel RKHS. Each iteration loops over k folds.
#
# Intermediate results are checkpointed after every iteration so that
# interrupted runs are recoverable.
#
# Usage:
#   Rscript RKHS_CV1.R <trait> [marker_file]
#   e.g.  Rscript RKHS_CV1.R DTF_blue
#   e.g.  Rscript RKHS_CV1.R DTF_blue auspak_for_rkhs.raw
#
# Output:
#   cv_results_CV1_<trait>_RKHS.csv   (final, all iterations)
#
# Seed scheme: 1000 + iteration (identical to BayesC for comparability).
#
# Reference:
#   Pérez & de los Campos (2014) Genetics 198:483-495
################################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript RKHS_CV1.R <trait> [marker_file]")
}

TRAIT       <- args[1]
MARKER_FILE <- if (length(args) >= 2) args[2] else "auspak_for_rkhs.raw"

cat("================================================================\n")
cat("RKHS CV1 — New Genotypes in Known Location-Years (Multi-Kernel)\n")
cat("Trait:", TRAIT, "\n")
cat("Marker file:", MARKER_FILE, "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

source("RKHS_utils.R")

# ── Load and prepare data ────────────────────────────────────────────────────
dat <- load_and_prepare_data(TRAIT, MARKER_FILE)
pheno_scaled <- dat$pheno
K_geno_list  <- dat$K_geno_list

# ── CV parameters ────────────────────────────────────────────────────────────
k_folds       <- 5
n_iterations  <- as.integer(Sys.getenv("RKHS_CV_ITERS", unset = "15"))
min_genotypes <- 10

# MCMC settings — match BayesC; override via environment variables for testing
nIter  <- as.integer(Sys.getenv("RKHS_NITER",  unset = "15000"))
burnIn <- as.integer(Sys.getenv("RKHS_BURNIN", unset = "5000"))
thin   <- as.integer(Sys.getenv("RKHS_THIN",   unset = "5"))

# ── Run CV1 ──────────────────────────────────────────────────────────────────

pheno_observed <- pheno_scaled  # unmasked copy
genotypes      <- unique(pheno_scaled$sample.id)
n_genotypes    <- length(genotypes)

cat("Data:", n_genotypes, "genotypes,",
    length(unique(pheno_scaled$location_year)), "location-years,",
    nrow(pheno_scaled), "observations\n")
cat("Settings:", k_folds, "folds,", n_iterations, "iterations\n\n")

checkpoint_file <- paste0("cv_results_CV1_", TRAIT, "_RKHS_checkpoint.csv")

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

# Resume from checkpoint if available
start_iter <- 1
if (file.exists(checkpoint_file)) {
  cv_results <- read.csv(checkpoint_file, stringsAsFactors = FALSE)
  completed_iters <- unique(cv_results$iteration)
  start_iter <- max(completed_iters) + 1
  cat("Resuming from checkpoint: iterations", paste(completed_iters, collapse = ","),
      "already done. Starting at iteration", start_iter, "\n\n")
}

if (start_iter <= n_iterations) {
  for (iter in start_iter:n_iterations) {
    current_seed <- 1000 + iter
    set.seed(current_seed)
    cat("Iteration", iter, "of", n_iterations, "(seed:", current_seed, ")\n")

    fold_assignments <- sample(rep(1:k_folds, length.out = n_genotypes))
    names(fold_assignments) <- genotypes

    for (fold in 1:k_folds) {
      cat("  Fold", fold, "of", k_folds, "\n")

      test_genotypes <- names(fold_assignments)[fold_assignments == fold]

      tryCatch({
        model_data <- pheno_scaled
        mask_rows  <- which(model_data$sample.id %in% test_genotypes)
        model_data[[TRAIT]][mask_rows] <- NA

        if (sum(!is.na(model_data[[TRAIT]])) < 50) {
          warning("Too few training observations"); next
        }

        saveAt <- paste0("RKHS_CV1_", TRAIT, "_i", iter, "_f", fold, "_")

        pred_values <- fit_rkhs_and_predict(
          model_data  = model_data,
          trait       = TRAIT,
          K_geno_list = K_geno_list,
          nIter = nIter, burnIn = burnIn, thin = thin, saveAt = saveAt
        )

        if (is.null(pred_values)) next

        eval_out <- evaluate_per_location_year(
          observed_data = pheno_observed,
          pred_values   = pred_values,
          test_indices  = mask_rows,
          trait = TRAIT, min_genotypes = min_genotypes
        )

        ly_results <- eval_out$metrics
        if (nrow(ly_results) > 0) {
          ly_results$iteration <- iter
          ly_results$fold      <- fold
          ly_results$seed      <- current_seed
          ly_results$cv_scheme <- "CV1"
          cv_results <- rbind(cv_results, ly_results)
        }

        if (nrow(eval_out$predictions) > 0) {
          preds <- eval_out$predictions
          preds$trait     <- TRAIT
          preds$iteration <- iter
          preds$fold      <- fold
          preds$seed      <- current_seed
          preds$cv_scheme <- "CV1"
          cv_predictions <- rbind(cv_predictions, preds)
        }

      }, error = function(e) {
        warning(paste("Error CV1 iter", iter, "fold", fold, ":", e$message))
      })
    }

    # Checkpoint after each iteration
    write.csv(cv_results, checkpoint_file, row.names = FALSE)
    cat("  Checkpoint saved after iteration", iter,
        "(", nrow(cv_results), "evaluations so far)\n\n")
  }
} else {
  cat("All iterations already completed.\n")
}

# ── Save final results ───────────────────────────────────────────────────────

outfile <- paste0("cv_results_CV1_", TRAIT, "_RKHS.csv")
write.csv(cv_results, outfile, row.names = FALSE)

pred_outfile <- paste0("predictions_CV1_", TRAIT, "_RKHS.csv")
write.csv(cv_predictions, pred_outfile, row.names = FALSE)

# Clean up checkpoint
if (file.exists(checkpoint_file)) file.remove(checkpoint_file)

cat("\n================================================================\n")
cat("CV1 complete for trait:", TRAIT, "\n")
cat("Total evaluations:", nrow(cv_results), "\n")
cat("Results saved to:", outfile, "\n")
cat("Predictions saved to:", pred_outfile, "(", nrow(cv_predictions), "rows )\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n")