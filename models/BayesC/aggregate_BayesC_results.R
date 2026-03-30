################################################################################
# aggregate_BayesC_results.R
#
# Collects per-iteration CV result CSVs into combined output files.
# Run manually after all SLURM jobs have completed.
#
# Usage:
#   Rscript aggregate_BayesC_results.R <trait>
#   Rscript aggregate_BayesC_results.R all       # process all traits
#
# Expects files in the working directory matching these patterns:
#   cv_results_CV1_<trait>_iter<NN>_BayesC.csv    (one per CV1 iteration)
#   cv_results_CV2_<trait>_iter<NN>_BayesC.csv    (one per CV2 iteration)
#   cv_results_CV0_<trait>_BayesC.csv             (one file)
#   cv_results_CrossLoc_<trait>_BayesC.csv        (one file)
#   predictions_CV1_<trait>_iter<NN>_BayesC.csv   (one per CV1 iteration)
#   predictions_CV2_<trait>_iter<NN>_BayesC.csv   (one per CV2 iteration)
#   predictions_CV0_<trait>_BayesC.csv            (one file)
#   predictions_CrossLoc_<trait>_BayesC.csv       (one file)
#
# Output:
#   cv_results_<trait>_BayesC_all_schemes.csv     (combined row-level results)
#   cv_summary_<trait>_BayesC_all_schemes.csv     (summary statistics)
#   predictions_<trait>_BayesC_all_schemes.csv    (combined individual predictions)
################################################################################

library(dplyr)

# ── Helper: read and combine CSVs matching a pattern ─────────────────────────

collect_csvs <- function(pattern, label) {
  files <- sort(list.files(pattern = pattern))
  if (length(files) == 0) {
    cat("  WARNING: No files found for", label, "(pattern:", pattern, ")\n")
    return(NULL)
  }
  cat("  ", label, ":", length(files), "files found\n")
  for (f in files) cat("    ", f, "\n")
  
  dfs <- lapply(files, function(f) {
    tryCatch(read.csv(f, stringsAsFactors = FALSE),
             error = function(e) {
               cat("    ERROR reading", f, ":", e$message, "\n")
               return(NULL)
             })
  })
  dfs <- dfs[!sapply(dfs, is.null)]
  
  if (length(dfs) == 0) return(NULL)
  
  combined <- bind_rows(dfs)
  cat("    -> Combined:", nrow(combined), "rows\n")
  return(combined)
}

# ── Helper: compute summary statistics (matches GBLUP pipeline) ─────────────

compute_summary <- function(cv_results) {
  if (is.null(cv_results) || nrow(cv_results) == 0) return(NULL)
  
  # For CV1/CV2: summarise by trait × location (across iterations/folds)
  # For CV0: summarise by trait × location (across location-years)
  # For CrossLoc: summarise by trait × cv_scheme (direction)
  
  summary_stats <- cv_results %>%
    group_by(trait, cv_scheme, location) %>%
    summarise(
      across(c(pearson, spearman, ndcg_at_10),
             list(mean = ~mean(., na.rm = TRUE),
                  sd   = ~sd(., na.rm = TRUE),
                  min  = ~min(., na.rm = TRUE),
                  max  = ~max(., na.rm = TRUE)),
             .names = "{.col}_{.fn}"),
      n_evaluations = n(),
      mean_n_test_genotypes = mean(n_test_genotypes, na.rm = TRUE),
      .groups = 'drop'
    )
  
  return(summary_stats)
}

# ── Main execution (skipped when sourced for testing) ─────────────────────────

if (sys.nframe() == 0) {

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript aggregate_BayesC_results.R <trait|all>")
}

VALID_TRAITS <- c('DTF_blue', 'DTH_blue', 'PtHt_blue', 'PcleLng_blue',
                  'SdLen_blue', 'TGW_blue', 'SdW_z_blue')

if (args[1] == "all") {
  traits_to_process <- VALID_TRAITS
} else {
  if (!args[1] %in% VALID_TRAITS) {
    stop("Invalid trait: '", args[1], "'\n",
         "  Valid: ", paste(VALID_TRAITS, collapse = ", "), ", or 'all'")
  }
  traits_to_process <- args[1]
}

# ── Process each trait ───────────────────────────────────────────────────────

for (TRAIT in traits_to_process) {
  
  cat("\n================================================================\n")
  cat("Aggregating results for trait:", TRAIT, "\n")
  cat("================================================================\n\n")
  
  all_parts <- list()
  any_missing <- FALSE
  
  # CV1: per-iteration files
  cv1_pattern <- paste0("^cv_results_CV1_", TRAIT, "_iter[0-9]+_BayesC\\.csv$")
  cv1 <- collect_csvs(cv1_pattern, "CV1")
  if (!is.null(cv1)) all_parts[["CV1"]] <- cv1
  else { cat("  -> CV1 MISSING\n"); any_missing <- TRUE }
  
  # CV2: per-iteration files
  cv2_pattern <- paste0("^cv_results_CV2_", TRAIT, "_iter[0-9]+_BayesC\\.csv$")
  cv2 <- collect_csvs(cv2_pattern, "CV2")
  if (!is.null(cv2)) all_parts[["CV2"]] <- cv2
  else { cat("  -> CV2 MISSING\n"); any_missing <- TRUE }
  
  # CV0: single file
  cv0_pattern <- paste0("^cv_results_CV0_", TRAIT, "_BayesC\\.csv$")
  cv0 <- collect_csvs(cv0_pattern, "CV0")
  if (!is.null(cv0)) all_parts[["CV0"]] <- cv0
  else { cat("  -> CV0 MISSING\n"); any_missing <- TRUE }
  
  # CrossLoc: single file
  xl_pattern <- paste0("^cv_results_CrossLoc_", TRAIT, "_BayesC\\.csv$")
  xl <- collect_csvs(xl_pattern, "CrossLoc")
  if (!is.null(xl)) all_parts[["CrossLoc"]] <- xl
  else { cat("  -> CrossLoc MISSING\n"); any_missing <- TRUE }
  
  if (any_missing) {
    cat("\n  WARNING: Some CV schemes have no results for", TRAIT, "\n")
    cat("  Aggregating what is available...\n")
  }
  
  if (length(all_parts) == 0) {
    cat("  ERROR: No results found for", TRAIT, "— skipping.\n")
    next
  }
  
  # Combine all schemes
  all_results <- bind_rows(all_parts)
  
  cat("\n  Combined results:", nrow(all_results), "rows\n")
  cat("  Breakdown:\n")
  print(table(all_results$cv_scheme))
  
  # Verify expected iteration counts
  if (!is.null(cv1)) {
    n_cv1_iters <- length(unique(cv1$iteration))
    cat("\n  CV1: found", n_cv1_iters, "iterations\n")
    if (n_cv1_iters < 15) {
      cat("  WARNING: Expected 15 CV1 iterations, found", n_cv1_iters, "\n")
      cat("  Missing iterations:",
          paste(setdiff(1:15, unique(cv1$iteration)), collapse = ", "), "\n")
    }
  }
  if (!is.null(cv2)) {
    n_cv2_iters <- length(unique(cv2$iteration))
    cat("  CV2: found", n_cv2_iters, "iterations\n")
    if (n_cv2_iters < 15) {
      cat("  WARNING: Expected 15 CV2 iterations, found", n_cv2_iters, "\n")
      cat("  Missing iterations:",
          paste(setdiff(1:15, unique(cv2$iteration)), collapse = ", "), "\n")
    }
  }
  
  # Save combined results
  outfile_results <- paste0("cv_results_", TRAIT, "_BayesC_all_schemes.csv")
  write.csv(all_results, outfile_results, row.names = FALSE)
  cat("\n  Results saved to:", outfile_results, "\n")
  
  # Compute and save summary
  summary_stats <- compute_summary(all_results)
  if (!is.null(summary_stats)) {
    outfile_summary <- paste0("cv_summary_", TRAIT, "_BayesC_all_schemes.csv")
    write.csv(summary_stats, outfile_summary, row.names = FALSE)
    cat("  Summary saved to:", outfile_summary, "\n")
    
    cat("\n  === Summary ===\n")
    print(as.data.frame(summary_stats[, c('trait', 'cv_scheme', 'location',
                                           'pearson_mean', 'pearson_sd',
                                           'spearman_mean', 'spearman_sd',
                                           'ndcg_at_10_mean', 'ndcg_at_10_sd',
                                           'n_evaluations')]))
  }
  
  # ── Aggregate individual predictions ──────────────────────────────────────
  cat("\n  --- Individual predictions ---\n")
  pred_parts <- list()

  pred_cv1_pat <- paste0("^predictions_CV1_", TRAIT, "_iter[0-9]+_BayesC\\.csv$")
  pred_cv1 <- collect_csvs(pred_cv1_pat, "CV1 predictions")
  if (!is.null(pred_cv1)) pred_parts[["CV1"]] <- pred_cv1

  pred_cv2_pat <- paste0("^predictions_CV2_", TRAIT, "_iter[0-9]+_BayesC\\.csv$")
  pred_cv2 <- collect_csvs(pred_cv2_pat, "CV2 predictions")
  if (!is.null(pred_cv2)) pred_parts[["CV2"]] <- pred_cv2

  pred_cv0_pat <- paste0("^predictions_CV0_", TRAIT, "_BayesC\\.csv$")
  pred_cv0 <- collect_csvs(pred_cv0_pat, "CV0 predictions")
  if (!is.null(pred_cv0)) pred_parts[["CV0"]] <- pred_cv0

  pred_xl_pat <- paste0("^predictions_CrossLoc_", TRAIT, "_BayesC\\.csv$")
  pred_xl <- collect_csvs(pred_xl_pat, "CrossLoc predictions")
  if (!is.null(pred_xl)) pred_parts[["CrossLoc"]] <- pred_xl

  if (length(pred_parts) > 0) {
    all_predictions <- bind_rows(pred_parts)
    outfile_preds <- paste0("predictions_", TRAIT, "_BayesC_all_schemes.csv")
    write.csv(all_predictions, outfile_preds, row.names = FALSE)
    cat("\n  Predictions saved to:", outfile_preds, "(", nrow(all_predictions), "rows )\n")
  } else {
    cat("\n  No prediction files found for", TRAIT, "\n")
  }

  # Optionally clean up per-iteration files
  cat("\n  Per-iteration files are preserved. To clean up, run:\n")
  cat("    rm cv_results_CV1_", TRAIT, "_iter*_BayesC.csv\n", sep = "")
  cat("    rm cv_results_CV2_", TRAIT, "_iter*_BayesC.csv\n", sep = "")
  cat("    rm predictions_CV1_", TRAIT, "_iter*_BayesC.csv\n", sep = "")
  cat("    rm predictions_CV2_", TRAIT, "_iter*_BayesC.csv\n", sep = "")
}

cat("\n================================================================\n")
cat("Aggregation complete.\n")
cat("Processed traits:", paste(traits_to_process, collapse = ", "), "\n")
cat("================================================================\n")

} # end if (sys.nframe() == 0)