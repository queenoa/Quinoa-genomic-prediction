# GBLUP cross-validation with z-score transformation
# Extended with multiple CV schemes:
#   CV1: New genotypes in known location-years
#   CV2: Sparse testing (known genotypes, incomplete location-years)
#   CV0: Leave-one-location-year-out
#   Cross-location transferability: within-location GBLUP -> predict other location
#
# Heterogeneous variance version:
#   residual = ~ dsum(~ units | location)
#   random   = ~ vm(sample.id, Ginv_sparse) + at(location):year
# Variance components from summary(model)$varcomp are captured per fit.
#
# REVISION NOTES:
# - fit_gblup_and_predict uses a unified NA-masking approach:
#   the full dataset is passed to asreml with target cells set to NA,
#   so all factor levels (locations, years) are represented in the model.
#   Used by CV0, CV1, and CV2.
# - Cross-location uses a separate, simpler within-location GBLUP
#   (intercept + G matrix only) because location-level fixed effects
#   cannot be estimated when all data for one location is NA. Within a
#   single location dsum/at(location) collapse to homogeneous, so the
#   cross-location model retains ~ units residual and ~ vm() random.

source("GBLUP_utils.R")


# ============================================================================
# CV1: New genotypes in known location-years
# Test genotypes have their trait values set to NA across all location-years.
# ============================================================================

run_cv1 <- function(pheno_data, Ginv_sparse, traits,
                     k_folds = 5, n_iterations = 15,
                     apply_zscore = TRUE, min_genotypes = 10) {

  cat("\n==========================================\n")
  cat("CV1: New genotypes in known location-years\n")
  cat("==========================================\n\n")
  cat("Preparing data...\n")

  pheno_data <- align_genotypes_to_gmatrix(pheno_data, Ginv_sparse)

  cat("Data: ", length(unique(pheno_data$sample.id)), "genotypes, ",
      length(levels(pheno_data$location_year)), "location-years, ",
      nrow(pheno_data), "observations\n")
  cat("Z-score:", ifelse(apply_zscore, "ON", "OFF"), "\n")
  cat("Settings:", k_folds, "folds,", n_iterations, "iterations,",
      "min genotypes per location-year:", min_genotypes, "\n\n")

  if (apply_zscore) {
    for (trait in traits) {
      pheno_data <- apply_location_year_scaling(pheno_data, trait)
    }
  }

  # Store unmasked data for evaluation
  pheno_data_observed <- pheno_data

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

  cv_varcomps <- data.frame()

  genotypes <- levels(pheno_data$sample.id)
  n_genotypes <- length(genotypes)

  cat("Starting CV:", n_genotypes, "genotypes |", k_folds, "folds |", n_iterations, "iterations\n\n")

  for (iter in 1:n_iterations) {
    current_seed <- 1000 + iter
    set.seed(current_seed)

    cat("Iteration", iter, "of", n_iterations, "(seed:", current_seed, ")\n")

    fold_assignments <- sample(rep(1:k_folds, length.out = n_genotypes))
    names(fold_assignments) <- genotypes

    for (fold in 1:k_folds) {
      cat("  Fold", fold, "of", k_folds, "\n")

      test_genotypes <- names(fold_assignments)[fold_assignments == fold]

      for (trait in traits) {
        cat("    Trait:", trait, "\n")

        tryCatch({

          # Mask test genotypes
          model_data <- pheno_data
          mask_rows <- model_data$sample.id %in% test_genotypes
          model_data[[trait]][mask_rows] <- NA

          model_out <- fit_gblup_and_predict(
            model_data = model_data,
            trait = trait,
            Ginv_sparse = Ginv_sparse
          )

          if (is.null(model_out)) next

          pred_values <- model_out$pred_values

          if (!is.null(model_out$varcomp) && nrow(model_out$varcomp) > 0) {
            vc <- model_out$varcomp
            vc$trait <- trait
            vc$iteration <- iter
            vc$fold <- fold
            vc$seed <- current_seed
            vc$cv_scheme <- "CV1"
            cv_varcomps <- rbind(cv_varcomps, vc)
          }

          eval_out <- evaluate_per_location_year(
            observed_data = pheno_data_observed,
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
            ly_results$cv_scheme <- "CV1"
            cv_results <- rbind(cv_results, ly_results)
          }

          if (nrow(eval_out$predictions) > 0) {
            preds <- eval_out$predictions
            preds$trait <- trait
            preds$iteration <- iter
            preds$fold <- fold
            preds$seed <- current_seed
            preds$cv_scheme <- "CV1"
            cv_predictions <- rbind(cv_predictions, preds)
          }

        }, error = function(e) {
          warning(paste("Error in trait", trait, ":", e$message))
        })
      }
    }
  }

  summary_stats <- summarise_cv_results(cv_results, "CV1")

  return(list(
    results = cv_results,
    predictions = cv_predictions,
    varcomps = cv_varcomps,
    summary = summary_stats,
    n_genotypes = n_genotypes,
    n_location_years = length(levels(pheno_data$location_year)),
    zscore_applied = apply_zscore,
    cv_scheme = "CV1"
  ))
}


# ============================================================================
# CV2: Sparse testing (5-fold stratified by location-year)
# Observed cells are assigned to 5 folds within each location-year so that
# each fold has balanced representation across environments. Each fold is
# held out in turn. Every observation is tested exactly once per iteration.
# Matches BayesC/RKHS CV2 design.
#
# NOTE: Cannot use evaluate_per_location_year() here because CV2 masks at
# the cell level — the same genotype may be masked in one location-year
# but observed in another. The join-by-genotype-ID approach would leak
# unmasked observations. Instead we track masked cells explicitly.
# ============================================================================

run_cv2 <- function(pheno_data, Ginv_sparse, traits,
                    k_folds = 5, n_iterations = 15,
                    apply_zscore = TRUE, min_genotypes = 10) {

  cat("\n==================================================\n")
  cat("CV2: Sparse testing (5-fold stratified)\n")
  cat("==================================================\n\n")
  cat("Preparing data...\n")

  pheno_data <- align_genotypes_to_gmatrix(pheno_data, Ginv_sparse)

  cat("Data: ", length(unique(pheno_data$sample.id)), "genotypes, ",
      length(levels(pheno_data$location_year)), "location-years, ",
      nrow(pheno_data), "observations\n")
  cat("Z-score:", ifelse(apply_zscore, "ON", "OFF"), "\n")
  cat("Settings:", k_folds, "folds,", n_iterations, "iterations,",
      "min genotypes per location-year:", min_genotypes, "\n\n")

  if (apply_zscore) {
    for (trait in traits) {
      pheno_data <- apply_location_year_scaling(pheno_data, trait)
    }
  }

  # Store unmasked data for evaluation
  pheno_data_observed <- pheno_data

  # Location lookup for adding location to predictions
  ly_loc <- pheno_data %>%
    mutate(location_year = as.character(location_year),
           location = as.character(location)) %>%
    distinct(location_year, location)

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

  cv_varcomps <- data.frame()

  for (iter in 1:n_iterations) {
    current_seed <- 2000 + iter
    set.seed(current_seed)
    cat("Iteration", iter, "of", n_iterations, "(seed:", current_seed, ")\n")

    for (trait in traits) {
      cat("  Trait:", trait, "\n")

      tryCatch({

        observed_idx <- which(!is.na(pheno_data[[trait]]))
        n_observed <- length(observed_idx)

        if (n_observed < 100) {
          warning(paste("Too few observations for", trait, ":", n_observed))
          next
        }

        # Stratified fold assignment: within each location-year, assign
        # observed cells to folds so each fold has balanced representation
        fold_assignments <- integer(nrow(pheno_data))  # 0 for unobserved
        location_years_obs <- as.character(pheno_data$location_year[observed_idx])

        for (ly in unique(location_years_obs)) {
          ly_obs_idx <- observed_idx[location_years_obs == ly]
          n_ly <- length(ly_obs_idx)
          fold_assignments[ly_obs_idx] <- sample(rep(1:k_folds, length.out = n_ly))
        }

        for (fold in 1:k_folds) {
          cat("    Fold", fold, "of", k_folds, "\n")

          mask_rows <- which(fold_assignments == fold)

          # Store observed values before masking
          masked_cells <- data.frame(
            sample.id = as.character(pheno_data$sample.id[mask_rows]),
            location_year = as.character(pheno_data$location_year[mask_rows]),
            observed = pheno_data_observed[[trait]][mask_rows],
            stringsAsFactors = FALSE
          )

          # Mask and fit
          model_data <- pheno_data
          model_data[[trait]][mask_rows] <- NA

          n_training <- sum(!is.na(model_data[[trait]]))
          if (n_training < 50) {
            warning("Too few training observations: ", n_training)
            next
          }

          # Safety check: report location-years with thin training data
          ly_train_counts <- table(
            as.character(pheno_data$location_year[which(!is.na(model_data[[trait]]))])
          )
          thin_lys <- names(ly_train_counts)[ly_train_counts < 10]
          if (length(thin_lys) > 0) {
            cat("      Warning: location-years with <10 training obs:",
                paste(thin_lys, collapse = ", "), "\n")
          }

          model_out <- fit_gblup_and_predict(
            model_data = model_data,
            trait = trait,
            Ginv_sparse = Ginv_sparse
          )

          if (is.null(model_out)) next

          pred_values <- model_out$pred_values

          if (!is.null(model_out$varcomp) && nrow(model_out$varcomp) > 0) {
            vc <- model_out$varcomp
            vc$trait <- trait
            vc$iteration <- iter
            vc$fold <- fold
            vc$seed <- current_seed
            vc$cv_scheme <- "CV2"
            cv_varcomps <- rbind(cv_varcomps, vc)
          }

          # Match predictions to the specifically masked cells
          pred_subset <- pred_values %>%
            select(sample.id, location_year, predicted = predicted.value)

          merged <- masked_cells %>%
            inner_join(pred_subset, by = c("sample.id", "location_year")) %>%
            filter(!is.na(predicted))

          if (nrow(merged) == 0) next

          # Add location info
          merged <- merged %>% left_join(ly_loc, by = "location_year")

          # Collect predictions
          if (nrow(merged) > 0) {
            preds <- data.frame(
              sample.id = merged$sample.id,
              location_year = merged$location_year,
              location = merged$location,
              observed = merged$observed,
              predicted = merged$predicted,
              trait = trait,
              iteration = iter,
              fold = fold,
              seed = current_seed,
              cv_scheme = "CV2",
              stringsAsFactors = FALSE
            )
            cv_predictions <- rbind(cv_predictions, preds)
          }

          # Evaluate per location-year
          for (ly in unique(merged$location_year)) {
            ly_data <- merged %>% filter(location_year == ly)
            n_geno <- length(unique(ly_data$sample.id))

            if (n_geno < min_genotypes) {
              cat("        Skipping", ly, "- only", n_geno, "test genotypes",
                  "(minimum:", min_genotypes, ")\n")
              next
            }

            eval_res <- evaluate_predictions(ly_data$observed, ly_data$predicted,
                                              trait_name = trait)

            if (!is.na(eval_res$pearson)) {
              loc <- ly_data$location[1]
              cv_results <- rbind(cv_results, data.frame(
                iteration = iter,
                fold = fold,
                trait = trait,
                location_year = ly,
                location = loc,
                pearson = eval_res$pearson,
                spearman = eval_res$spearman,
                ndcg_at_10 = eval_res$ndcg_at_10,
                seed = current_seed,
                n_test_genotypes = n_geno,
                cv_scheme = "CV2",
                stringsAsFactors = FALSE
              ))

              cat("        ", ly, "- r:", round(eval_res$pearson, 3),
                  "| rho:", round(eval_res$spearman, 3),
                  "| NDCG@10:", round(eval_res$ndcg_at_10, 3),
                  "| N:", n_geno, "\n")
            }
          }
        }

      }, error = function(e) {
        warning(paste("Error in trait", trait, ":", e$message))
      })
    }
  }

  summary_stats <- summarise_cv_results(cv_results, "CV2")

  return(list(
    results = cv_results,
    predictions = cv_predictions,
    varcomps = cv_varcomps,
    summary = summary_stats,
    zscore_applied = apply_zscore,
    cv_scheme = "CV2"
  ))
}


# ============================================================================
# CV0: Leave-one-location-year-out
# All rows kept; trait set to NA for the held-out location-year.
# ============================================================================

run_cv0 <- function(pheno_data, Ginv_sparse, traits,
                    apply_zscore = TRUE, min_genotypes = 10) {

  cat("\n=============================================\n")
  cat("CV0: Leave-one-location-year-out\n")
  cat("=============================================\n\n")
  cat("Preparing data...\n")

  pheno_data <- align_genotypes_to_gmatrix(pheno_data, Ginv_sparse)

  location_years <- levels(pheno_data$location_year)

  cat("Data: ", length(unique(pheno_data$sample.id)), "genotypes, ",
      length(location_years), "location-years, ",
      nrow(pheno_data), "observations\n")
  cat("Z-score:", ifelse(apply_zscore, "ON", "OFF"), "\n")
  cat("Min genotypes:", min_genotypes, "\n\n")

  if (apply_zscore) {
    for (trait in traits) {
      pheno_data <- apply_location_year_scaling(pheno_data, trait)
    }
  }

  pheno_data_observed <- pheno_data

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
    trait = character(), cv_scheme = character(),
    stringsAsFactors = FALSE
  )

  cv_varcomps <- data.frame()

  for (held_out_ly in location_years) {
    cat("Holding out:", held_out_ly, "\n")

    held_out_location <- pheno_data %>%
      filter(location_year == held_out_ly) %>%
      pull(location) %>%
      unique() %>%
      as.character()

    for (trait in traits) {
      cat("  Trait:", trait, "\n")

      tryCatch({

        test_genos_in_ly <- pheno_data_observed %>%
          filter(location_year == held_out_ly & !is.na(.data[[trait]])) %>%
          pull(sample.id) %>%
          unique() %>%
          as.character()

        if (length(test_genos_in_ly) < min_genotypes) {
          cat("    Skipping - only", length(test_genos_in_ly),
              "genotypes observed in held-out location-year\n")
          next
        }

        # Mask held-out location-year
        model_data <- pheno_data
        mask_rows <- model_data$location_year == held_out_ly
        model_data[[trait]][mask_rows] <- NA

        # CV0 uses homogeneous residuals: AUS traits with only two of three
        # years cause the heterogeneous (dsum / at(location):year) fit to
        # fail when a whole location-year is held out.
        model <- tryCatch({
          asreml(
            fixed = as.formula(paste(trait, "~ location")),
            random = ~ vm(sample.id, Ginv_sparse) + location:year,
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
          cat("    Model did not converge\n")
          next
        }

        predictions <- tryCatch({
          predict(model, classify = "sample.id:location:year",
                  pworkspace = 256e06, trace = FALSE)
        }, error = function(e) {
          warning(paste("Prediction error for", trait, ":", e$message))
          return(NULL)
        })

        if (is.null(predictions$pvals)) next

        pred_values <- predictions$pvals
        pred_values$sample.id <- as.character(pred_values$sample.id)
        pred_values$location <- as.character(pred_values$location)
        pred_values$year <- as.character(pred_values$year)
        pred_values$location_year <- paste(pred_values$location,
                                            pred_values$year, sep = "_")

        # Filter to only location-years present in the input data
        real_lys <- unique(as.character(model_data$location_year))
        pred_values <- pred_values[pred_values$location_year %in% real_lys, ]

        vc_raw <- summary(model)$varcomp
        if (!is.null(vc_raw) && nrow(vc_raw) > 0) {
          vc <- data.frame(
            component_name = rownames(vc_raw),
            vc_raw,
            row.names = NULL,
            stringsAsFactors = FALSE,
            check.names = FALSE
          )
          vc$trait <- trait
          vc$held_out_location_year <- held_out_ly
          vc$held_out_location <- held_out_location
          vc$cv_scheme <- "CV0"
          cv_varcomps <- rbind(cv_varcomps, vc)
        }

        # Extract predictions for the held-out location-year
        pred_in_ly <- pred_values %>%
          filter(location_year == held_out_ly & sample.id %in% test_genos_in_ly) %>%
          select(sample.id, location_year, predicted = predicted.value)

        obs_in_ly <- pheno_data_observed %>%
          mutate(sample.id = as.character(sample.id),
                 location_year = as.character(location_year)) %>%
          filter(location_year == held_out_ly &
                   sample.id %in% test_genos_in_ly &
                   !is.na(.data[[trait]])) %>%
          select(sample.id, location_year, observed = all_of(trait))

        merged <- obs_in_ly %>%
          inner_join(pred_in_ly, by = c("sample.id", "location_year")) %>%
          filter(!is.na(predicted))

        n_geno <- length(unique(merged$sample.id))

        # Collect predictions
        if (nrow(merged) > 0) {
          preds <- data.frame(
            sample.id = merged$sample.id,
            location_year = merged$location_year,
            location = held_out_location,
            observed = merged$observed,
            predicted = merged$predicted,
            trait = trait,
            cv_scheme = "CV0",
            stringsAsFactors = FALSE
          )
          cv_predictions <- rbind(cv_predictions, preds)
        }

        if (n_geno < min_genotypes) {
          cat("    Skipping - only", n_geno, "genotypes with predictions\n")
          next
        }

        eval_res <- evaluate_predictions(merged$observed, merged$predicted,
                                          trait_name = trait)

        if (!is.na(eval_res$pearson)) {
          cv_results <- rbind(cv_results, data.frame(
            iteration = NA,
            fold = NA,
            trait = trait,
            location_year = held_out_ly,
            location = held_out_location,
            pearson = eval_res$pearson,
            spearman = eval_res$spearman,
            ndcg_at_10 = eval_res$ndcg_at_10,
            seed = NA,
            n_test_genotypes = n_geno,
            cv_scheme = "CV0",
            stringsAsFactors = FALSE
          ))

          cat("    ", held_out_ly, "- r:", round(eval_res$pearson, 3),
              "| rho:", round(eval_res$spearman, 3),
              "| NDCG@10:", round(eval_res$ndcg_at_10, 3),
              "| N:", n_geno, "\n")
        }

      }, error = function(e) {
        warning(paste("Error in trait", trait, "for", held_out_ly, ":", e$message))
      })
    }
  }

  summary_stats <- summarise_cv_results(cv_results, "CV0")

  return(list(
    results = cv_results,
    predictions = cv_predictions,
    varcomps = cv_varcomps,
    summary = summary_stats,
    zscore_applied = apply_zscore,
    cv_scheme = "CV0"
  ))
}


# ============================================================================
# Cross-location transferability (simple within-location GBLUP)
#
# Fits trait ~ 1 + vm(sample.id, G) on the SOURCE location only.
# Extracts genomic BLUPs and correlates them with observed phenotypes
# in each TARGET location-year.
#
# No location or year fixed/random effects -- prediction into the other
# location comes purely from genomic relationships.
# ============================================================================

run_cross_location <- function(pheno_data, Ginv_sparse, traits,
                                apply_zscore = TRUE, min_genotypes = 10) {

  cat("\n=============================================\n")
  cat("Cross-location transferability\n")
  cat("(within-location GBLUP -> predict other location)\n")
  cat("=============================================\n\n")

  pheno_data <- align_genotypes_to_gmatrix(pheno_data, Ginv_sparse)

  locations <- unique(as.character(pheno_data$location))

  cat("Data:", length(unique(pheno_data$sample.id)), "genotypes,",
      length(levels(pheno_data$location_year)), "location-years,",
      nrow(pheno_data), "observations\n")
  cat("Locations:", paste(locations, collapse = ", "), "\n")
  cat("Z-score:", ifelse(apply_zscore, "ON", "OFF"), "\n")
  cat("Min genotypes:", min_genotypes, "\n\n")

  if (apply_zscore) {
    for (trait in traits) {
      pheno_data <- apply_location_year_scaling(pheno_data, trait)
    }
  }

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
    trait = character(), train_location = character(), cv_scheme = character(),
    stringsAsFactors = FALSE
  )

  cv_varcomps <- data.frame()

  for (train_loc in locations) {
    predict_loc <- setdiff(locations, train_loc)

    cat("Training on:", train_loc, "-> Predicting:",
        paste(predict_loc, collapse = ", "), "\n")

    train_data <- pheno_data %>% filter(location == train_loc)

    for (trait in traits) {
      cat("  Trait:", trait, "\n")

      tryCatch({

        n_train_obs <- sum(!is.na(train_data[[trait]]))
        if (n_train_obs < 50) {
          cat("    Skipping - only", n_train_obs, "training observations\n")
          next
        }

        # Fit simple within-location GBLUP (intercept + G matrix only)
        model <- asreml(
          fixed  = as.formula(paste(trait, "~ 1")),
          random = ~ vm(sample.id, Ginv_sparse),
          residual = ~ units,
          na.action = na.method(y = "include"),
          workspace = 256e06,
          data = train_data,
          maxit = 200,
          trace = FALSE
        )

        if (is.null(model) || !model$converge) {
          cat("    Model did not converge\n")
          next
        }

        # Capture variance components from this within-location fit
        vc_raw <- summary(model)$varcomp
        if (!is.null(vc_raw) && nrow(vc_raw) > 0) {
          vc <- data.frame(
            component_name = rownames(vc_raw),
            vc_raw,
            row.names = NULL,
            stringsAsFactors = FALSE,
            check.names = FALSE
          )
          vc$trait <- trait
          vc$train_location <- train_loc
          vc$cv_scheme <- paste0("CrossLoc_", train_loc, "->",
                                 paste(predict_loc, collapse = "+"))
          cv_varcomps <- rbind(cv_varcomps, vc)
        }

        # Extract genomic BLUPs
        blups <- predict(model, classify = "sample.id",
                         pworkspace = 256e06, trace = FALSE)

        if (is.null(blups$pvals)) {
          cat("    Could not extract BLUPs\n")
          next
        }

        gebv <- blups$pvals %>%
          mutate(sample.id = as.character(sample.id)) %>%
          select(sample.id, gebv = predicted.value)

        cat("    Extracted GEBVs for", nrow(gebv), "genotypes\n")

        # Evaluate in each target location-year
        target_lys <- pheno_data %>%
          filter(location %in% predict_loc) %>%
          pull(location_year) %>%
          unique() %>%
          as.character()

        cv_label <- paste0("CrossLoc_", train_loc, "->",
                           paste(predict_loc, collapse = "+"))

        for (ly in target_lys) {

          obs_in_ly <- pheno_data %>%
            filter(location_year == ly & !is.na(.data[[trait]])) %>%
            mutate(sample.id = as.character(sample.id),
                   location = as.character(location)) %>%
            select(sample.id, location_year, location,
                   observed = all_of(trait))

          merged <- obs_in_ly %>%
            inner_join(gebv, by = "sample.id") %>%
            filter(!is.na(gebv))

          n_geno <- length(unique(merged$sample.id))

          # Collect predictions
          if (nrow(merged) > 0) {
            preds <- data.frame(
              sample.id = merged$sample.id,
              location_year = as.character(merged$location_year),
              location = as.character(merged$location),
              observed = merged$observed,
              predicted = merged$gebv,
              trait = trait,
              train_location = train_loc,
              cv_scheme = cv_label,
              stringsAsFactors = FALSE
            )
            cv_predictions <- rbind(cv_predictions, preds)
          }

          if (n_geno < min_genotypes) {
            cat("      Skipping", ly, "- only", n_geno,
                "overlapping genotypes (min:", min_genotypes, ")\n")
            next
          }

          eval_res <- evaluate_predictions(
            merged$observed, merged$gebv, trait_name = trait
          )

          if (!is.na(eval_res$pearson)) {
            cv_results <- rbind(cv_results, data.frame(
              iteration = NA,
              fold = NA,
              trait = trait,
              location_year = ly,
              location = as.character(merged$location[1]),
              pearson = eval_res$pearson,
              spearman = eval_res$spearman,
              ndcg_at_10 = eval_res$ndcg_at_10,
              seed = NA,
              n_test_genotypes = n_geno,
              cv_scheme = cv_label,
              stringsAsFactors = FALSE
            ))

            cat("      ", ly, "- r:", round(eval_res$pearson, 3),
                "| rho:", round(eval_res$spearman, 3),
                "| NDCG@10:", round(eval_res$ndcg_at_10, 3),
                "| N:", n_geno, "\n")
          }
        }

      }, error = function(e) {
        warning(paste("Error in trait", trait, ":", e$message))
      })
    }
  }

  # Summary (group by direction, not location)
  if (nrow(cv_results) > 0) {
    summary_stats <- cv_results %>%
      group_by(trait, cv_scheme) %>%
      summarise(
        across(c(pearson, spearman, ndcg_at_10),
               list(mean = ~mean(., na.rm = TRUE),
                    sd   = ~sd(., na.rm = TRUE),
                    min  = ~min(., na.rm = TRUE),
                    max  = ~max(., na.rm = TRUE)),
               .names = "{.col}_{.fn}"),
        n_target_location_years = n(),
        mean_n_test_genotypes = mean(n_test_genotypes, na.rm = TRUE),
        .groups = 'drop'
      )

    cat("\n=== Cross-Location Transferability Summary ===\n")
    print(summary_stats[, c('trait', 'cv_scheme', 'pearson_mean', 'pearson_sd',
                             'spearman_mean', 'spearman_sd',
                             'ndcg_at_10_mean', 'ndcg_at_10_sd',
                             'n_target_location_years')])
  } else {
    warning("No successful cross-location results")
    summary_stats <- data.frame()
  }

  return(list(
    results = cv_results,
    predictions = cv_predictions,
    varcomps = cv_varcomps,
    summary = summary_stats,
    zscore_applied = apply_zscore,
    cv_scheme = "CrossLocation"
  ))
}


# ============================================================================
# WRAPPER: Run all CV schemes and combine results
# ============================================================================

run_all_cv_schemes <- function(pheno_data, Ginv_sparse, traits,
                                k_folds = 5, n_iterations = 15,
                                apply_zscore = TRUE, min_genotypes = 10) {

  cat("==========================================================\n")
  cat("Running all cross-validation schemes for GBLUP\n")
  cat("==========================================================\n\n")

  results <- list()

  results$cv1 <- run_cv1(pheno_data, Ginv_sparse, traits,
                          k_folds = k_folds, n_iterations = n_iterations,
                          apply_zscore = apply_zscore,
                          min_genotypes = min_genotypes)

  results$cv2 <- run_cv2(pheno_data, Ginv_sparse, traits,
                          k_folds = k_folds,
                          n_iterations = n_iterations,
                          apply_zscore = apply_zscore,
                          min_genotypes = min_genotypes)

  results$cv0 <- run_cv0(pheno_data, Ginv_sparse, traits,
                          apply_zscore = apply_zscore,
                          min_genotypes = min_genotypes)

  results$cross_loc <- run_cross_location(pheno_data, Ginv_sparse, traits,
                                           apply_zscore = apply_zscore,
                                           min_genotypes = min_genotypes)

  all_results <- bind_rows(
    results$cv1$results,
    results$cv2$results,
    results$cv0$results,
    results$cross_loc$results
  )

  cat("\n==========================================================\n")
  cat("All CV schemes complete\n")
  cat("Total evaluations:", nrow(all_results), "\n")
  cat("Breakdown:\n")
  print(table(all_results$cv_scheme))
  cat("==========================================================\n")

  results$all_results <- all_results
  return(results)
}


# ============================================================================
# RUN PIPELINE
# ============================================================================
# sys.nframe() == 0 means the script was invoked directly (Rscript / source
# from R prompt). When it's source()'d from inside a function (e.g. tests),
# the run-pipeline block is skipped so the script can be loaded for its
# function definitions only.

if (sys.nframe() == 0) {

# Load G-inverse (produced by G_matrix_GBLUP.R)
load("Ginv_sparse_GBLUP.RData")

# Load phenotype data
pheno_data <- read.csv(file.path("..", "..", "data",
                                  "AUSPAK_phenotypes_GP_input.csv"),
                        stringsAsFactors = FALSE)

# Run all CV schemes
all_cv <- run_all_cv_schemes(
  pheno_data = pheno_data,
  Ginv_sparse= Ginv_sparse,
  traits = c("DTF", "DTH", "PtHt", "PcleLng",
                  "SdLen", "TGW", "SdW_z"),
  k_folds = 5,
  n_iterations = 1,
  apply_zscore = TRUE,
  min_genotypes = 10
)


# ============================================================================
# SAVE OUTPUTS
# ============================================================================

scheme_labels <- c(cv1 = "CV1", cv2 = "CV2", cv0 = "CV0",
                   cross_loc = "CrossLoc")

traits <- c("DTF", "DTH", "PtHt", "PcleLng",
                  "SdLen", "TGW", "SdW_z")

# --- Per-scheme, per-trait CSVs ---
for (scheme in names(scheme_labels)) {
  label <- scheme_labels[[scheme]]
  res    <- all_cv[[scheme]]$results
  preds  <- all_cv[[scheme]]$predictions
  vcomps <- all_cv[[scheme]]$varcomps

  for (trait in traits) {
    trait_res <- res[res$trait == trait, ]
    if (nrow(trait_res) > 0) {
      outfile <- paste0("cv_results_", label, "_", trait, "_GBLUP.csv")
      write.csv(trait_res, outfile, row.names = FALSE)
      cat("Saved:", outfile, "(", nrow(trait_res), "rows)\n")
    }

    trait_preds <- preds[preds$trait == trait, ]
    if (nrow(trait_preds) > 0) {
      outfile <- paste0("predictions_", label, "_", trait, "_GBLUP.csv")
      write.csv(trait_preds, outfile, row.names = FALSE)
      cat("Saved:", outfile, "(", nrow(trait_preds), "rows)\n")
    }

    if (!is.null(vcomps) && nrow(vcomps) > 0) {
      trait_vcomps <- vcomps[vcomps$trait == trait, ]
      if (nrow(trait_vcomps) > 0) {
        outfile <- paste0("varcomps_", label, "_", trait, "_GBLUP.csv")
        write.csv(trait_vcomps, outfile, row.names = FALSE)
        cat("Saved:", outfile, "(", nrow(trait_vcomps), "rows)\n")
      }
    }
  }
}

# --- Combined all-schemes per-trait CSVs ---
all_preds <- bind_rows(
  all_cv$cv1$predictions,
  all_cv$cv2$predictions,
  all_cv$cv0$predictions,
  all_cv$cross_loc$predictions
)

all_summaries <- bind_rows(
  lapply(names(scheme_labels), function(scheme) {
    s <- all_cv[[scheme]]$summary
    if (nrow(s) > 0) { s$cv_scheme <- scheme_labels[[scheme]] }
    s
  })
)

all_varcomps <- bind_rows(
  all_cv$cv1$varcomps,
  all_cv$cv2$varcomps,
  all_cv$cv0$varcomps,
  all_cv$cross_loc$varcomps
)

for (trait in traits) {
  # All results across schemes
  trait_all <- all_cv$all_results[all_cv$all_results$trait == trait, ]
  if (nrow(trait_all) > 0) {
    outfile <- paste0("cv_results_", trait, "_GBLUP_all_schemes.csv")
    write.csv(trait_all, outfile, row.names = FALSE)
    cat("Saved:", outfile, "(", nrow(trait_all), "rows)\n")
  }

  # Summary statistics across schemes
  trait_summary <- all_summaries[all_summaries$trait == trait, ]
  if (nrow(trait_summary) > 0) {
    outfile <- paste0("cv_summary_", trait, "_GBLUP_all_schemes.csv")
    write.csv(trait_summary, outfile, row.names = FALSE)
    cat("Saved:", outfile, "(", nrow(trait_summary), "rows)\n")
  }

  # All predictions across schemes
  trait_preds <- all_preds[all_preds$trait == trait, ]
  if (nrow(trait_preds) > 0) {
    outfile <- paste0("predictions_", trait, "_GBLUP_all_schemes.csv")
    write.csv(trait_preds, outfile, row.names = FALSE)
    cat("Saved:", outfile, "(", nrow(trait_preds), "rows)\n")
  }

  # All variance components across schemes
  if (!is.null(all_varcomps) && nrow(all_varcomps) > 0) {
    trait_vcomps <- all_varcomps[all_varcomps$trait == trait, ]
    if (nrow(trait_vcomps) > 0) {
      outfile <- paste0("varcomps_", trait, "_GBLUP_all_schemes.csv")
      write.csv(trait_vcomps, outfile, row.names = FALSE)
      cat("Saved:", outfile, "(", nrow(trait_vcomps), "rows)\n")
    }
  }
}

cat("\n=== All outputs saved ===\n")

} # end if (sys.nframe() == 0)
