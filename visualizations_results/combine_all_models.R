################################################################################
# combine_all_models.R
#
# Combines cv_results_*_all_schemes.csv, cv_summary_*_all_schemes.csv, and
# predictions_*_all_schemes.csv from each model subdirectory into single
# files with a 'model' column.
#
# Usage (from results/):
#   Rscript combine_all_models.R
#
# Output:
#   cv_results_all_models.csv
#   cv_summary_all_models.csv
#   predictions_all_models.csv.gz   (gzipped — ~5-10x smaller than raw CSV;
#                                    pandas reads it transparently with
#                                    pd.read_csv(..., compression='infer'))
################################################################################

library(dplyr)

# Model directories and their display names
model_dirs <- c(
  "GBLUP_heterogeneous"          = "GBLUP",
  "BayesC"         = "BayesC",
  "RKHS"           = "RKHS",
  "LightGBM"       = "LightGBM"
)

# ── Combine cv_results ───────────────────────────────────────────────────────

cat("=== Combining cv_results files ===\n\n")

all_results <- list()

for (dir_name in names(model_dirs)) {
  model_label <- model_dirs[[dir_name]]
  files <- list.files(dir_name,
                      pattern = "^cv_results_.*_all_schemes\\.csv$",
                      full.names = TRUE)

  if (length(files) == 0) {
    cat("  WARNING: No cv_results files in", dir_name, "\n")
    next
  }

  cat("  ", model_label, ":", length(files), "files\n")

  dfs <- lapply(files, function(f) {
    df <- read.csv(f, stringsAsFactors = FALSE)
    df$model <- model_label
    df
  })

  all_results[[dir_name]] <- bind_rows(dfs)
}

combined_results <- bind_rows(all_results)
cat("\nTotal rows:", nrow(combined_results), "\n")
cat("Models:", paste(unique(combined_results$model), collapse = ", "), "\n")
cat("Traits:", paste(sort(unique(combined_results$trait)), collapse = ", "), "\n")
cat("CV schemes:", paste(sort(unique(combined_results$cv_scheme)), collapse = ", "), "\n")

write.csv(combined_results, "cv_results_all_models.csv", row.names = FALSE)
cat("\nSaved: cv_results_all_models.csv\n")

# ── Combine cv_summary ───────────────────────────────────────────────────────

cat("\n=== Combining cv_summary files ===\n\n")

all_summaries <- list()

for (dir_name in names(model_dirs)) {
  model_label <- model_dirs[[dir_name]]
  files <- list.files(dir_name,
                      pattern = "^cv_summary_.*_all_schemes\\.csv$",
                      full.names = TRUE)

  if (length(files) == 0) {
    cat("  WARNING: No cv_summary files in", dir_name, "\n")
    next
  }

  cat("  ", model_label, ":", length(files), "files\n")

  dfs <- lapply(files, function(f) {
    df <- read.csv(f, stringsAsFactors = FALSE)
    df$model <- model_label
    df
  })

  all_summaries[[dir_name]] <- bind_rows(dfs)
}

combined_summary <- bind_rows(all_summaries)
write.csv(combined_summary, "cv_summary_all_models.csv", row.names = FALSE)
cat("\nSaved: cv_summary_all_models.csv\n")

# ── Combine predictions ──────────────────────────────────────────────────────
# Predictions are larger than results/summaries (one row per test observation),
# so we write gzipped CSV. CV0/CrossLoc files lack iteration/fold/seed columns
# and CrossLoc has an extra train_location column — bind_rows fills missing
# columns with NA so the union schema is preserved.

cat("\n=== Combining predictions files ===\n\n")

all_predictions <- list()

for (dir_name in names(model_dirs)) {
  model_label <- model_dirs[[dir_name]]
  files <- list.files(dir_name,
                      pattern = "^predictions_.*_all_schemes\\.csv$",
                      full.names = TRUE)

  if (length(files) == 0) {
    cat("  WARNING: No predictions files in", dir_name, "\n")
    next
  }

  cat("  ", model_label, ":", length(files), "files\n")

  dfs <- lapply(files, function(f) {
    df <- read.csv(f, stringsAsFactors = FALSE)
    df$model <- model_label
    df
  })

  all_predictions[[dir_name]] <- bind_rows(dfs)
}

combined_predictions <- bind_rows(all_predictions)
cat("\nTotal rows:", format(nrow(combined_predictions), big.mark = ","), "\n")
cat("Columns:", paste(names(combined_predictions), collapse = ", "), "\n")
cat("Models:", paste(unique(combined_predictions$model), collapse = ", "), "\n")
cat("CV schemes:",
    paste(sort(unique(combined_predictions$cv_scheme)), collapse = ", "), "\n")

pred_path <- "predictions_all_models.csv.gz"
gz <- gzfile(pred_path, "w")
write.csv(combined_predictions, gz, row.names = FALSE)
close(gz)
cat("\nSaved:", pred_path,
    sprintf("(%.1f MB)\n", file.info(pred_path)$size / 1024^2))

cat("\n=== Done ===\n")
