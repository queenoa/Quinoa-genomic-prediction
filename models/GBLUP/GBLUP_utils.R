################################################################################
# GBLUP_utils.R — Shared utility functions for GBLUP multi-environment CV
#
# Sourced by:
#   GBLUP.R
#
# Contains:
#   - Global constants (traits, lower-is-better)
#   - NDCG@k calculation
#   - Evaluation metrics (Pearson, Spearman, NDCG@10)
#   - Z-score standardisation by location-year
#   - Genotype alignment with G matrix
#   - GBLUP model fitting and prediction (ASReml)
#   - Per-location-year evaluation
#   - CV results summarisation
#
# Functions shared with BayesC/RKHS (identical logic):
#   VALID_TRAITS, LOWER_IS_BETTER_TRAITS
#   apply_location_year_scaling()
#   calculate_ndcg()
#   evaluate_predictions()
#   summarise_cv_results()
#
# Functions specific to GBLUP:
#   align_genotypes_to_gmatrix()  — factor-level alignment for ASReml G matrix
#   fit_gblup_and_predict()       — ASReml GBLUP fitting + classify prediction
#   evaluate_per_location_year()  — join-based (test_ids) rather than index-based
################################################################################

library(dplyr)
library(asreml)

# ── Global constants ─────────────────────────────────────────────────────────

VALID_TRAITS <- c('DTF', 'DTH', 'PtHt', 'PcleLng',
                  'SdLen', 'TGW', 'SdW_z')

LOWER_IS_BETTER_TRAITS <- c('DTF', 'DTH', 'PtHt')

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

# ── Validate and align genotype factor levels with G matrix ──────────────────

align_genotypes_to_gmatrix <- function(pheno_data, Ginv_sparse) {

  pheno_data$sample.id <- as.factor(as.character(pheno_data$sample.id))
  pheno_data$location_year <- as.factor(as.character(pheno_data$location_year))
  pheno_data$location <- as.factor(as.character(pheno_data$location))
  pheno_data$year <- as.factor(as.character(pheno_data$year))

  G_rownames <- rownames(Ginv_sparse)
  if (is.null(G_rownames)) G_rownames <- attr(Ginv_sparse, "rowNames")

  geno_in_data <- unique(as.character(pheno_data$sample.id))
  missing_from_G <- geno_in_data[!geno_in_data %in% G_rownames]
  if (length(missing_from_G) > 0) {
    stop(paste("These genotypes are in phenotype data but not in G matrix:",
               paste(head(missing_from_G, 10), collapse = ", ")))
  }

  pheno_data$sample.id <- factor(pheno_data$sample.id, levels = G_rownames)

  if (!identical(levels(pheno_data$sample.id), G_rownames)) {
    stop("Failed to align genotype factor levels with G matrix order")
  }

  cat("Genotype alignment verified:", length(geno_in_data),
      "genotypes in data,", length(G_rownames), "in G matrix\n")

  return(pheno_data)
}

# ── Fit GBLUP model and return predictions ───────────────────────────────────
#
# Takes the full dataset with target cells already set to NA.
# All factor levels remain in the model so ASReml can predict for any
# location-year combination, including held-out ones.
# Used by CV0, CV1, and CV2 (not cross-location).
#
# Heterogeneous variance structure:
#   residual = ~ dsum(~ units | location)   — separate residual var per location
#   random   = ~ vm(sample.id, Ginv_sparse) + at(location):year
#                                          — separate year-within-location var
#                                            for each location
#
# dsum requires data sorted by the sectioning factor (location), so we sort
# inside this function before fitting.
#
# Returns a list:
#   pred_values — data.frame of predicted values per sample.id × location × year
#   varcomp     — data.frame of variance components from summary(model)$varcomp,
#                 with the rownames as a `component_name` column

fit_gblup_and_predict <- function(model_data, trait, Ginv_sparse) {

  non_na_count <- sum(!is.na(model_data[[trait]]))
  if (non_na_count < 50) {
    warning(paste("Insufficient training data for", trait, ":", non_na_count, "obs"))
    return(NULL)
  }

  # dsum() in the residual term requires data sorted by the sectioning factor
  model_data <- model_data[order(as.character(model_data$location)), ]

  model <- tryCatch({
    asreml(
      fixed = as.formula(paste(trait, "~ location")),
      random = ~ vm(sample.id, Ginv_sparse) + at(location):year,
      residual = ~ dsum(~ units | location),
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

  predictions <- tryCatch({
    predict(model, classify = "sample.id:location:year",
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

  # Coerce factor columns to character to avoid factor-level mismatch
  # in downstream joins. ASReml predict with classify produces the full
  # factorial (e.g., PAK_2017, AUS_2020) which must be filtered out.
  pred_values$sample.id <- as.character(pred_values$sample.id)
  pred_values$location <- as.character(pred_values$location)
  pred_values$year <- as.character(pred_values$year)
  pred_values$location_year <- paste(pred_values$location, pred_values$year, sep = "_")

  # Filter to only location-years that exist in the input data
  real_lys <- unique(as.character(model_data$location_year))
  pred_values <- pred_values[pred_values$location_year %in% real_lys, ]

  vc <- summary(model)$varcomp
  varcomp_df <- data.frame(
    component_name = rownames(vc),
    vc,
    row.names = NULL,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  return(list(
    pred_values = pred_values,
    varcomp = varcomp_df
  ))
}

# ── Evaluate per location-year ───────────────────────────────────────────────
#
# Join-based approach: takes test_ids (character vector of genotype names)
# rather than integer row indices (as in BayesC/RKHS). The join produces the
# same result but matches the ASReml predict output format.
#
# Returns list(metrics, predictions) matching BayesC/RKHS interface.

evaluate_per_location_year <- function(observed_data, pred_values, test_ids,
                                        trait, min_genotypes = 10) {

  results <- data.frame(
    trait = character(), location_year = character(), location = character(),
    pearson = numeric(), spearman = numeric(), ndcg_at_10 = numeric(),
    n_test_genotypes = integer(), stringsAsFactors = FALSE
  )

  test_ids <- as.character(test_ids)

  obs <- observed_data %>%
    mutate(sample.id = as.character(sample.id),
           location_year = as.character(location_year),
           location = as.character(location)) %>%
    filter(sample.id %in% test_ids & !is.na(.data[[trait]])) %>%
    select(sample.id, location_year, location, observed = all_of(trait))

  pred <- pred_values %>%
    filter(sample.id %in% test_ids) %>%
    select(sample.id, location_year, predicted = predicted.value)

  merged <- obs %>%
    inner_join(pred, by = c("sample.id", "location_year")) %>%
    filter(!is.na(predicted))

  if (nrow(merged) == 0) return(list(metrics = results, predictions = merged))

  for (ly in unique(merged$location_year)) {
    ly_data <- merged %>% filter(location_year == ly)
    n_geno <- length(unique(ly_data$sample.id))

    if (n_geno < min_genotypes) {
      cat("      Skipping", ly, "- only", n_geno, "test genotypes",
          "(minimum:", min_genotypes, ")\n")
      next
    }

    eval_res <- evaluate_predictions(ly_data$observed, ly_data$predicted,
                                      trait_name = trait)

    if (!is.na(eval_res$pearson)) {
      loc <- ly_data$location[1]
      results <- rbind(results, data.frame(
        trait = trait, location_year = ly, location = loc,
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

# ── Summarise CV results ─────────────────────────────────────────────────────

summarise_cv_results <- function(cv_results, scheme_name) {
  if (nrow(cv_results) == 0) {
    warning(paste("No successful cross-validation results for", scheme_name))
    return(data.frame())
  }

  summary_stats <- cv_results %>%
    group_by(trait, location) %>%
    summarise(
      across(c(pearson, spearman, ndcg_at_10),
             list(mean = ~mean(., na.rm = TRUE),
                  sd = ~sd(., na.rm = TRUE),
                  min = ~min(., na.rm = TRUE),
                  max = ~max(., na.rm = TRUE)),
             .names = "{.col}_{.fn}"),
      n_evaluations = n(),
      mean_n_test_genotypes = mean(n_test_genotypes, na.rm = TRUE),
      .groups = 'drop'
    )

  cat("\n===", scheme_name, "Cross-Validation Summary ===\n")
  print(summary_stats[, c('trait', 'location', 'pearson_mean', 'pearson_sd',
                           'spearman_mean', 'spearman_sd',
                           'ndcg_at_10_mean', 'ndcg_at_10_sd', 'n_evaluations')])

  return(summary_stats)
}

cat("GBLUP_utils.R loaded successfully.\n")
