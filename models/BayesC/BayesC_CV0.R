################################################################################
# BayesC_CV0.R
#
# Runs CV0 (leave-one-location-year-out) for one trait.
# Deterministic single-pass scheme — no random folds.
# Results and predictions are saved incrementally after each location-year.
#
# Usage:
#   Rscript BayesC_CV0.R <trait> [marker_file]
#   e.g.  Rscript BayesC_CV0.R DTF_blue
#
# Output:
#   cv_results_CV0_<trait>_BayesC.csv
#   predictions_CV0_<trait>_BayesC.csv
################################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript BayesC_CV0.R <trait> [marker_file]")
}

TRAIT       <- args[1]
MARKER_FILE <- if (length(args) >= 2) args[2] else "pruned05_AUSPAK_for_bayesC.raw"

cat("================================================================\n")
cat("BayesC CV0 — Leave-one-location-year-out\n")
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
outfile_results <- paste0("cv_results_CV0_", TRAIT, "_BayesC.csv")
outfile_preds   <- paste0("predictions_CV0_", TRAIT, "_BayesC.csv")

# ── CV0: Leave-one-location-year-out ────────────────────────────────────────

pheno_observed <- pheno_scaled
location_years <- sort(unique(pheno_scaled$location_year))

cat("Data:", length(unique(pheno_scaled$sample.id)), "genotypes,",
    length(location_years), "location-years,",
    nrow(pheno_scaled), "observations\n\n")

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
  trait = character(), cv_scheme = character(),
  stringsAsFactors = FALSE
)
write.csv(preds_header, outfile_preds, row.names = FALSE)

n_evaluations <- 0

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

    saveAt <- paste0("BayesC_CV0_", TRAIT, "_", held_out_ly, "_")

    pred_values <- fit_bayesc_and_predict(
      model_data = model_data, trait = TRAIT, X_geno = X_geno,
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
      write.table(preds_df, outfile_preds, append = TRUE, sep = ",",
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
      write.table(results_df, outfile_results, append = TRUE, sep = ",",
                  row.names = FALSE, col.names = FALSE, quote = TRUE)
      n_evaluations <- n_evaluations + 1

      cat("  ", held_out_ly, "- r:", round(eval_res$pearson, 3),
          "| rho:", round(eval_res$spearman, 3),
          "| NDCG@10:", round(eval_res$ndcg_at_10, 3),
          "| N:", n_geno, "\n")
    }

  }, error = function(e) {
    warning(paste("Error CV0 for", held_out_ly, ":", e$message))
  })
}

cat("\n================================================================\n")
cat("CV0 complete for trait:", TRAIT, "\n")
cat("CV0 evaluations:", n_evaluations, "\n")
cat("Results:", outfile_results, "\n")
cat("Predictions:", outfile_preds, "\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n")
