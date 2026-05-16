# GBLUP CV1 fitted SEPARATELY by location
#
# Same CV1 design as the global model in GBLUP.R (5-fold CV on genotypes,
# 15 iterations, seed = 1000+iter, z-score by location-year), but the model
# is fit independently on each location's subset.
#
# Goal: compare prediction accuracies of single-location models against the
# pooled global model that uses data from both locations.
#
# Model structure (per location):
#   fixed    = trait ~ 1                    (no location effect — single location)
#   random   = vm(sample.id, Ginv_sparse) + year
#   residual = ~ units
#
# Genotypes are sampled from those observed in the focal location only;
# predictions and evaluation are restricted to that location's location-years.
#
# Additionally fits a single-location GBLUP per (location, trait) on the full
# data (no CV masking) using the same model structure, extracts genetic and
# residual variance components plus h2, and writes them to
# per_location_variances_GBLUP.csv.

source("GBLUP_utils.R")


# ============================================================================
# Per-location GBLUP fit and prediction
# ============================================================================

fit_gblup_per_location <- function(model_data, trait, Ginv_sparse) {

  non_na_count <- sum(!is.na(model_data[[trait]]))
  if (non_na_count < 50) {
    warning(paste("Insufficient training data for", trait, ":", non_na_count, "obs"))
    return(NULL)
  }

  n_years <- length(unique(as.character(model_data$year)))

  random_formula <- if (n_years > 1) {
    ~ vm(sample.id, Ginv_sparse) + year
  } else {
    ~ vm(sample.id, Ginv_sparse)
  }

  model <- tryCatch({
    asreml(
      fixed    = as.formula(paste(trait, "~ 1")),
      random   = random_formula,
      residual = ~ units,
      na.action = na.method(y = "include"),
      workspace = 256e06,
      data = model_data,
      maxit = 200,
      trace = FALSE
    )
  }, error = function(e) {
    warning(paste("Model fitting error for", trait, ":", e$message))
    return(NULL)
  })

  if (is.null(model) || !model$converge) {
    warning(paste("Model convergence failed for", trait))
    return(NULL)
  }

  classify_str <- if (n_years > 1) "sample.id:year" else "sample.id"

  predictions <- tryCatch({
    predict(model, classify = classify_str,
            pworkspace = 256e06, trace = FALSE)
  }, error = function(e) {
    warning(paste("Prediction error for", trait, ":", e$message))
    return(NULL)
  })

  if (is.null(predictions$pvals)) {
    warning(paste("Could not extract predictions for trait", trait))
    return(NULL)
  }

  pred_values <- predictions$pvals
  pred_values$sample.id <- as.character(pred_values$sample.id)

  loc_label <- as.character(model_data$location[1])

  if (n_years > 1) {
    pred_values$year <- as.character(pred_values$year)
  } else {
    pred_values$year <- as.character(model_data$year[1])
  }
  pred_values$location <- loc_label
  pred_values$location_year <- paste(loc_label, pred_values$year, sep = "_")

  real_lys <- unique(as.character(model_data$location_year))
  pred_values <- pred_values[pred_values$location_year %in% real_lys, ]

  return(pred_values)
}


# ============================================================================
# Variance components: full-data per-location GBLUP (no CV masking)
# Same model structure as fit_gblup_per_location, but fit on the full
# location subset and returning (genetic, residual) variances and h2.
# ============================================================================

fit_single_location_gblup <- function(loc_data, trait, Ginv_sparse) {

  n_obs <- sum(!is.na(loc_data[[trait]]))
  if (n_obs < 30) {
    warning(paste("Insufficient data for", trait, ":", n_obs, "obs"))
    return(NULL)
  }

  n_years <- length(unique(as.character(loc_data$year)))

  random_formula <- if (n_years > 1) {
    ~ vm(sample.id, Ginv_sparse) + year
  } else {
    ~ vm(sample.id, Ginv_sparse)
  }

  model <- tryCatch({
    asreml(
      fixed    = as.formula(paste(trait, "~ 1")),
      random   = random_formula,
      residual = ~ units,
      na.action = na.method(y = "include"),
      workspace = 256e06,
      data = loc_data,
      maxit = 200,
      trace = FALSE
    )
  }, error = function(e) {
    warning(paste("Model fitting error for", trait, ":", e$message))
    return(NULL)
  })

  if (is.null(model) || !model$converge) {
    warning(paste("Model convergence failed for", trait))
    return(NULL)
  }

  vc <- summary(model)$varcomp

  res_row <- grep("units!R|units!units", rownames(vc), value = TRUE)
  if (length(res_row) == 0) res_row <- rownames(vc)[nrow(vc)]

  gen_row <- grep("vm\\(sample.id", rownames(vc), value = TRUE)

  res_var <- vc[res_row[1], "component"]
  gen_var <- if (length(gen_row) > 0) vc[gen_row[1], "component"] else NA
  h2 <- if (!is.na(gen_var)) gen_var / (gen_var + res_var) else NA

  return(data.frame(
    n_obs = n_obs,
    genetic_variance = gen_var,
    residual_variance = res_var,
    h2 = h2,
    stringsAsFactors = FALSE
  ))
}


run_per_location_variances <- function(pheno_data, Ginv_sparse, traits,
                                        apply_zscore = TRUE) {

  cat("\n==========================================================\n")
  cat("Per-location GBLUP variance component estimation\n")
  cat("==========================================================\n\n")

  pheno_data <- align_genotypes_to_gmatrix(pheno_data, Ginv_sparse)

  if (apply_zscore) {
    for (trait in traits) {
      pheno_data <- apply_location_year_scaling(pheno_data, trait)
    }
  }

  locations <- levels(factor(pheno_data$location))

  cat("Data:", length(unique(pheno_data$sample.id)), "genotypes,",
      length(locations), "locations,",
      nrow(pheno_data), "observations\n")
  cat("Z-score:", ifelse(apply_zscore, "ON", "OFF"), "\n\n")

  variance_results <- data.frame(
    trait = character(),
    location = character(),
    n_obs = integer(),
    genetic_variance = numeric(),
    residual_variance = numeric(),
    h2 = numeric(),
    stringsAsFactors = FALSE
  )

  for (loc in locations) {
    cat("Location:", loc, "\n")

    loc_data <- pheno_data %>% filter(location == loc)

    if (nrow(loc_data) < 30) {
      cat("  Skipping - too few observations\n")
      next
    }

    for (trait in traits) {
      cat("  Trait:", trait, "")

      vc <- fit_single_location_gblup(loc_data, trait, Ginv_sparse)

      if (is.null(vc)) {
        cat(" -> failed\n")
        next
      }

      cat(" -> Vg:", round(vc$genetic_variance, 4),
          "| Ve:", round(vc$residual_variance, 4),
          "| h2:", round(vc$h2, 3),
          "| N:", vc$n_obs, "\n")

      variance_results <- rbind(variance_results, data.frame(
        trait = trait,
        location = loc,
        n_obs = vc$n_obs,
        genetic_variance = vc$genetic_variance,
        residual_variance = vc$residual_variance,
        h2 = vc$h2,
        stringsAsFactors = FALSE
      ))
    }
  }

  return(variance_results)
}


# ============================================================================
# CV1 per location: new genotypes in known location-years (single location)
# ============================================================================

run_cv1_per_location <- function(pheno_data, Ginv_sparse, traits,
                                  k_folds = 5, n_iterations = 15,
                                  apply_zscore = TRUE, min_genotypes = 10) {

  cat("\n==================================================\n")
  cat("CV1 per-location (separate GBLUP for each location)\n")
  cat("==================================================\n\n")

  pheno_data <- align_genotypes_to_gmatrix(pheno_data, Ginv_sparse)

  if (apply_zscore) {
    for (trait in traits) {
      pheno_data <- apply_location_year_scaling(pheno_data, trait)
    }
  }

  pheno_data_observed <- pheno_data

  locations <- sort(unique(as.character(pheno_data$location)))

  cv_results <- data.frame(
    iteration = integer(), fold = integer(), trait = character(),
    location_year = character(), location = character(),
    pearson = numeric(), spearman = numeric(), ndcg_at_10 = numeric(),
    seed = integer(), n_test_genotypes = integer(),
    cv_scheme = character(), train_location = character(),
    stringsAsFactors = FALSE
  )

  cv_predictions <- data.frame(
    sample.id = character(), location_year = character(),
    location = character(), observed = numeric(), predicted = numeric(),
    trait = character(), iteration = integer(), fold = integer(),
    seed = integer(), cv_scheme = character(), train_location = character(),
    stringsAsFactors = FALSE
  )

  for (loc in locations) {
    cat("\n========== Location:", loc, "==========\n\n")

    loc_data <- pheno_data %>% filter(location == loc)
    loc_data_obs <- pheno_data_observed %>% filter(location == loc)

    loc_genotypes <- unique(as.character(loc_data$sample.id))
    n_loc_genos <- length(loc_genotypes)

    cat("Data:", n_loc_genos, "genotypes,",
        length(unique(as.character(loc_data$location_year))), "location-years,",
        nrow(loc_data), "observations\n")
    cat("Z-score:", ifelse(apply_zscore, "ON", "OFF"), "\n")
    cat("Settings:", k_folds, "folds,", n_iterations, "iterations,",
        "min genotypes per location-year:", min_genotypes, "\n\n")

    for (iter in 1:n_iterations) {
      current_seed <- 1000 + iter
      set.seed(current_seed)

      cat("Iteration", iter, "of", n_iterations, "(seed:", current_seed, ")\n")

      fold_assignments <- sample(rep(1:k_folds, length.out = n_loc_genos))
      names(fold_assignments) <- loc_genotypes

      for (fold in 1:k_folds) {
        cat("  Fold", fold, "of", k_folds, "\n")

        test_genotypes <- names(fold_assignments)[fold_assignments == fold]

        for (trait in traits) {
          cat("    Trait:", trait, "\n")

          tryCatch({

            model_data <- loc_data
            mask_rows <- model_data$sample.id %in% test_genotypes
            model_data[[trait]][mask_rows] <- NA

            pred_values <- fit_gblup_per_location(
              model_data = model_data,
              trait = trait,
              Ginv_sparse = Ginv_sparse
            )

            if (is.null(pred_values)) next

            eval_out <- evaluate_per_location_year(
              observed_data = loc_data_obs,
              pred_values = pred_values,
              test_ids = test_genotypes,
              trait = trait,
              min_genotypes = min_genotypes
            )

            ly_results <- eval_out$metrics
            if (nrow(ly_results) > 0) {
              ly_results$iteration <- iter
              ly_results$fold <- fold
              ly_results$seed <- current_seed
              ly_results$cv_scheme <- paste0("CV1_", loc)
              ly_results$train_location <- loc
              cv_results <- rbind(cv_results, ly_results)
            }

            if (nrow(eval_out$predictions) > 0) {
              preds <- eval_out$predictions
              preds$trait <- trait
              preds$iteration <- iter
              preds$fold <- fold
              preds$seed <- current_seed
              preds$cv_scheme <- paste0("CV1_", loc)
              preds$train_location <- loc
              cv_predictions <- rbind(cv_predictions, preds)
            }

          }, error = function(e) {
            warning(paste("Error in trait", trait, "(", loc, "):", e$message))
          })
        }
      }
    }
  }

  summary_stats <- summarise_cv_results(cv_results, "CV1 per-location")

  return(list(
    results = cv_results,
    predictions = cv_predictions,
    summary = summary_stats,
    zscore_applied = apply_zscore,
    cv_scheme = "CV1_per_location"
  ))
}


# ============================================================================
# RUN PIPELINE
# ============================================================================
# sys.nframe() == 0 means the script was invoked directly (Rscript / source
# from R prompt). When it's source()'d from inside a function (e.g. tests),
# the run-pipeline block is skipped so the script can be loaded for its
# function definitions only.

if (sys.nframe() == 0) {

load("Ginv_sparse_GBLUP.RData")

pheno_data <- read.csv(file.path("..", "..", "data",
                                  "AUSPAK_phenotypes_GP_input.csv"),
                        stringsAsFactors = FALSE)

traits <- c("DTF", "DTH", "PtHt", "PcleLng",
            "SdLen", "TGW", "SdW_z")


# --- Per-location variance components (full-data fit, no CV) ---

pheno_for_vc <- pheno_data
pheno_for_vc$year <- as.factor(as.character(pheno_for_vc$year))
pheno_for_vc$location <- as.factor(as.character(pheno_for_vc$location))
pheno_for_vc$sample.id <- as.factor(as.character(pheno_for_vc$sample.id))

variance_results <- run_per_location_variances(
  pheno_data   = pheno_for_vc,
  Ginv_sparse  = Ginv_sparse,
  traits       = traits,
  apply_zscore = TRUE
)

write.csv(variance_results, "per_location_variances_GBLUP.csv", row.names = FALSE)
cat("\nSaved: per_location_variances_GBLUP.csv (",
    nrow(variance_results), "rows)\n")


# --- CV1 per location ---

cv1_per_loc <- run_cv1_per_location(
  pheno_data   = pheno_data,
  Ginv_sparse  = Ginv_sparse,
  traits       = traits,
  k_folds      = 5,
  n_iterations = 15,
  apply_zscore = TRUE,
  min_genotypes = 10
)


# ============================================================================
# SAVE OUTPUTS
# ============================================================================

# Per-trait, per-location CSVs
res   <- cv1_per_loc$results
preds <- cv1_per_loc$predictions
locations <- sort(unique(as.character(res$train_location)))

for (loc in locations) {
  for (trait in traits) {
    trait_res <- res[res$trait == trait & res$train_location == loc, ]
    if (nrow(trait_res) > 0) {
      outfile <- paste0("cv_results_CV1_", loc, "_", trait, "_GBLUP.csv")
      write.csv(trait_res, outfile, row.names = FALSE)
      cat("Saved:", outfile, "(", nrow(trait_res), "rows)\n")
    }

    trait_preds <- preds[preds$trait == trait & preds$train_location == loc, ]
    if (nrow(trait_preds) > 0) {
      outfile <- paste0("predictions_CV1_", loc, "_", trait, "_GBLUP.csv")
      write.csv(trait_preds, outfile, row.names = FALSE)
      cat("Saved:", outfile, "(", nrow(trait_preds), "rows)\n")
    }
  }
}

# Combined results across locations (per trait)
for (trait in traits) {
  trait_res <- res[res$trait == trait, ]
  if (nrow(trait_res) > 0) {
    outfile <- paste0("cv_results_", trait, "_GBLUP_CV1_per_location.csv")
    write.csv(trait_res, outfile, row.names = FALSE)
    cat("Saved:", outfile, "(", nrow(trait_res), "rows)\n")
  }

  trait_preds <- preds[preds$trait == trait, ]
  if (nrow(trait_preds) > 0) {
    outfile <- paste0("predictions_", trait, "_GBLUP_CV1_per_location.csv")
    write.csv(trait_preds, outfile, row.names = FALSE)
    cat("Saved:", outfile, "(", nrow(trait_preds), "rows)\n")
  }
}

# Summary across all traits and locations
if (nrow(cv1_per_loc$summary) > 0) {
  write.csv(cv1_per_loc$summary,
            "cv_summary_GBLUP_CV1_per_location.csv",
            row.names = FALSE)
  cat("Saved: cv_summary_GBLUP_CV1_per_location.csv\n")
}

cat("\n=== All outputs saved ===\n")

} # end if (sys.nframe() == 0)
