library(dplyr)
library(asreml)

# ---- Diagnostic helpers ---------------------------------------------------

check_model_status <- function(model) {
  converged <- model$converge
  stable_parameters <- TRUE

  if (!is.null(model$warn.list)) {
    for (warn in model$warn.list) {
      if (grepl("changed by more than", warn, fixed = TRUE)) {
        stable_parameters <- FALSE
        break
      }
    }
  }

  list(converged = converged, stable = stable_parameters)
}

# Q-Q correlation is used instead of Shapiro-Wilk because SW is over-sensitive
# on large n.
check_residual_normality <- function(residuals, qq_correlation_threshold = 0.95) {
  if (sum(!is.na(residuals)) < 3) {
    return(list(
      correlation = NA_real_,
      severe_deviation = TRUE,
      interpretation = "Insufficient non-NA residuals for normality check"
    ))
  }

  residuals <- residuals[!is.na(residuals)]

  tryCatch({
    qq_data <- qqnorm(residuals, plot = FALSE)
    qq_correlation <- cor(qq_data$x, qq_data$y, use = "complete.obs")

    if (is.na(qq_correlation)) {
      return(list(
        correlation = NA_real_,
        severe_deviation = TRUE,
        interpretation = "Unable to calculate Q-Q correlation"
      ))
    }

    interpretation <- if (qq_correlation >= 0.98) {
      "Excellent normality"
    } else if (qq_correlation >= 0.95) {
      "Acceptable normality"
    } else if (qq_correlation >= 0.90) {
      "Moderate deviation from normality"
    } else {
      "Substantial deviation from normality"
    }

    list(
      correlation = qq_correlation,
      severe_deviation = qq_correlation < qq_correlation_threshold,
      interpretation = interpretation
    )
  }, error = function(e) {
    list(
      correlation = NA_real_,
      severe_deviation = TRUE,
      interpretation = paste("Error in normality check:", e$message)
    )
  })
}

gen_heritability <- function(model) {
  vc_m2 <- summary(model)$varcomp
  hv <- which(row.names(vc_m2) == "accession")
  vv <- vc_m2[hv, 1]
  hh <- coefficients(model)$random
  hh <- data.frame(hh)
  hh$std.error <- sqrt(model$sigma2 * model$vcoeff$random)
  hh <- hh[grep("accession_*", dimnames(hh)[[1]]), ]
  1 - mean(hh[, "std.error"])^2 / vv
}
