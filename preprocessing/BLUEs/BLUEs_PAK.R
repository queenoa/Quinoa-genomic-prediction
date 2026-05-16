# One BLUE per accession per trait across all PAK years.
# Matches AUS (year as random); only year 2021-22 is replicated, so
# trial_replicate enters the model only at that level via at().
# Run from preprocessing/BLUEs/:  Rscript BLUEs_PAK.R

source("BLUEs_utils.R")

TRAITS <- c("DTF", "DTH", "PtHt", "PcleLng", "SdLen", "TGW", "SdW_z")
LOCATION <- "PAK"
OUTPUT_DIR <- "BLUEs_results/PAK"
REPLICATED_YEAR <- "2021-22"
MAX_ITERATIONS <- 20
QQ_CORRELATION_THRESHOLD <- 0.95

pak <- read.csv(
  "../../data/PAK_phenotypes_raw.csv",
  header = TRUE,
  na.strings = c("", "NA"),
  stringsAsFactors = FALSE
)
pak$year <- factor(pak$year)
pak$accession <- factor(pak$accession)
pak$trial_replicate <- factor(pak$trial_replicate)

stopifnot(REPLICATED_YEAR %in% levels(pak$year))

message(sprintf("[%s] rows: %d", LOCATION, nrow(pak)))
message(sprintf("[%s] rows per year:", LOCATION))
print(table(pak$year))

sequenced_accessions <- pak %>%
  filter(!is.na(SampleName) & SampleName != "") %>%
  select(accession, SampleName) %>%
  distinct()
message(sprintf("[%s] sequenced accessions (with SampleName): %d",
                LOCATION, nrow(sequenced_accessions)))

plot_dir <- file.path(OUTPUT_DIR, "diagnostic_plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

all_results <- list()
heritability_results <- data.frame(
  Trait = character(),
  Heritability = numeric(),
  Model_Converged = logical(),
  stringsAsFactors = FALSE
)

for (trait in TRAITS) {
  tryCatch({
    message(sprintf("[%s] processing trait: %s", LOCATION, trait))

    model_data <- pak %>%
      select(all_of(c("accession", "year", "trial_replicate", "SampleName", trait))) %>%
      filter(!is.na(!!sym(trait)))

    if (nrow(model_data) == 0) {
      message(sprintf("[%s] no data for %s — skipping", LOCATION, trait))
      next
    }

    names(model_data)[names(model_data) == trait] <- "response"

    model_blues <- asreml(
      fixed = response ~ accession,
      random = ~ year + at(year, "2021-22"):trial_replicate,
      residual = ~ units,
      data = model_data,
      na.action = na.method(y = "exclude", x = "exclude"),
      trace = FALSE,
      workspace = "1gb"
    )
    iter_blues <- 1
    status_blues <- check_model_status(model_blues)
    while (iter_blues <= MAX_ITERATIONS && (!status_blues$converged || !status_blues$stable)) {
      model_blues <- update(model_blues)
      status_blues <- check_model_status(model_blues)
      iter_blues <- iter_blues + 1
    }

    model_h2 <- asreml(
      fixed = response ~ 1,
      random = ~ accession + year + at(year, "2021-22"):trial_replicate,
      residual = ~ units,
      data = model_data,
      na.action = na.method(y = "exclude", x = "exclude"),
      trace = FALSE,
      workspace = "1gb"
    )
    iter_h2 <- 1
    status_h2 <- check_model_status(model_h2)
    while (iter_h2 <= MAX_ITERATIONS && (!status_h2$converged || !status_h2$stable)) {
      model_h2 <- update(model_h2)
      status_h2 <- check_model_status(model_h2)
      iter_h2 <- iter_h2 + 1
    }

    h2 <- tryCatch(
      gen_heritability(model_h2),
      error = function(e) {
        message(sprintf("[%s] heritability error for %s: %s",
                        LOCATION, trait, e$message))
        NA_real_
      }
    )

    heritability_results <- rbind(
      heritability_results,
      data.frame(
        Trait = trait,
        Heritability = h2,
        Model_Converged = status_h2$converged && status_h2$stable
      )
    )

    residuals <- resid(model_blues)
    normality_check <- check_residual_normality(residuals, QQ_CORRELATION_THRESHOLD)

    pred <- predict(model_blues,
                    classify = "accession",
                    vcov = TRUE,
                    sed = TRUE,
                    pworkspace = 64e6)
    pred_df <- pred$pvals

    if (is.null(pred_df) || nrow(pred_df) == 0) {
      message(sprintf("[%s] no predictions for %s — skipping", LOCATION, trait))
      next
    }

    pred_df_sequenced <- pred_df %>%
      filter(accession %in% sequenced_accessions$accession) %>%
      left_join(sequenced_accessions, by = "accession")

    message(sprintf("[%s] %s predictions (sequenced only): %d",
                    LOCATION, trait, nrow(pred_df_sequenced)))

    blues_df <- data.frame(
      trait = trait,
      location = LOCATION,
      accession = pred_df_sequenced$accession,
      SampleName = pred_df_sequenced$SampleName,
      BLUE = pred_df_sequenced$predicted.value,
      SE = pred_df_sequenced$std.error,
      Convergence_Status = ifelse(
        status_blues$converged && status_blues$stable,
        "Fully converged", "Check convergence"
      ),
      Iterations = iter_blues - 1,
      QQ_Correlation = normality_check$correlation,
      Normality_Status = normality_check$interpretation
    )

    pdf(file.path(plot_dir, paste0(trait, "_diagnostics.pdf")))
    par(mfrow = c(2, 2))
    plot(fitted(model_blues), residuals,
         main = paste0(trait, "\nResiduals vs Fitted"),
         xlab = "Fitted values", ylab = "Residuals")
    abline(h = 0, col = "red", lty = 2)
    qqnorm(residuals,
           main = paste0("Normal Q-Q Plot\nCorrelation: ",
                         round(normality_check$correlation, 3)))
    qqline(residuals, col = "red")
    hist(residuals,
         main = paste0("Histogram of residuals\n",
                       normality_check$interpretation),
         breaks = 30)
    boxplot(response ~ year, data = model_data,
            main = "Distribution by Year",
            xlab = "Year", ylab = trait)
    dev.off()

    all_results[[trait]] <- blues_df
  }, error = function(e) {
    message(sprintf("[%s] error in %s: %s", LOCATION, trait, e$message))
    all_results[[trait]] <- NULL
  })
}

if (length(all_results) > 0) {
  final_results <- do.call(rbind, all_results)
  write.csv(final_results,
            file = file.path(OUTPUT_DIR, "trait_BLUEs_results.csv"),
            row.names = FALSE)
} else {
  message(sprintf("[%s] no results generated", LOCATION))
}

write.csv(heritability_results,
          file = file.path(OUTPUT_DIR, "trait_heritability_results.csv"),
          row.names = FALSE)
