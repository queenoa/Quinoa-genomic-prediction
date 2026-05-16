# BLUEs pipeline

Computes one BLUE per accession per location per trait for use in phenotype
summaries (correlation plots, descriptive statistics). These BLUEs are **not**
the input to the genomic prediction models — those use raw per-location-year
means for every environment, including PAK 2021-22. PAK 2021-22 is the only
replicated trial, but without spatial correction or additional covariates the
BLUEs do not differ meaningfully from raw means, so means are used throughout
to keep the GP input uniform.

## Files

- `BLUEs_utils.R` — shared diagnostic helpers only (`check_model_status`,
  `check_residual_normality`, `gen_heritability`). Loads `dplyr` and `asreml`.
- `BLUEs_AUS.R` — self-contained: loads `../../data/AUS_phenotypes_raw.csv`
  and runs the full per-trait loop with AUS formulas inline.
- `BLUEs_PAK.R` — self-contained: loads `../../data/PAK_phenotypes_raw.csv`
  and runs the full per-trait loop with PAK formulas inline (including the
  `at(year, "2021-22"):trial_replicate` term).
- `merge_BLUEs_for_phenotype_plots.R` — concatenates
  `BLUEs_results/AUS/trait_BLUEs_results.csv` and
  `BLUEs_results/PAK/trait_BLUEs_results.csv` into
  `BLUEs_results/AUSPAK_BLUEs_for_phenotype_plots.csv` (long format, one row
  per trait × accession × location). This file is the input for correlation
  plots and descriptive phenotype summaries.
- `../build_GP_phenotype_input.R` (in `preprocessing/`, not this directory) —
  builds the phenotype input file for the genomic prediction models at
  `data/AUSPAK_phenotypes_GP_input.csv` using raw per-accession means for
  every location-year (AUS 2017/2018/2019, PAK 2019-20, PAK 2020-21,
  PAK 2021-22) computed directly from `data/{AUS,PAK}_phenotypes_raw.csv`.
  Wide format — one row per (location, year, accession), one column per
  trait. Only sequenced accessions (non-NA `SampleName`) are kept. The
  pipeline does not consume any BLUEs outputs.

The per-trait loop is duplicated across the two scripts by design.
`asreml()` does NSE on its `fixed` / `random` arguments, so formulas passed
through a wrapper function get rejected as `"object is of mode call, please
simplify"`. Keeping each script self-contained with literal formulas
sidesteps this entirely.

## Models

Both locations produce one BLUE per accession (collapsed across years within
a location):

| Location | BLUEs `fixed` | BLUEs `random` | H2 `fixed` | H2 `random` |
|----------|---------------|----------------|------------|-------------|
| AUS (all years) | `response ~ accession` | `~ year` | `response ~ 1` | `~ accession + year` |
| PAK (all years) | `response ~ accession` | `~ year + at(year, "2021-22"):trial_replicate` | `response ~ 1` | `~ accession + year + at(year, "2021-22"):trial_replicate` |

Predictions come from the BLUEs model via `predict(classify = "accession")`.
Heritability is computed from the companion random-accession model using
`gen_heritability()`.

### Why PAK uses `at(year, "2021-22"):trial_replicate`

Only the 2021-22 PAK trial is replicated (two reps). 2019-20 and 2020-21 have
a single rep each, so a global `year:trial_replicate` term would be confounded
with `year` for those years. `at(year, "2021-22")` restricts the
replicate term to the single year where it carries information, and keeps the
rest of the model identical to AUS.

## Fit all accessions, filter predictions at the end

Each script intentionally fits the ASReml model on the **full** phenotype
dataset — including accessions that were never sequenced (no `SampleName`).
Only the returned prediction table is filtered to sequenced accessions.

Why: `year` and residual variance components are estimated more stably when
every observed plot contributes. Dropping unsequenced accessions before
fitting throws away information that informs the random year effect and the
residual, without changing which accessions we actually want BLUEs for.
Filtering happens once at the end, on `pred$pvals`, via
`filter(accession %in% sequenced_accessions$accession)`.

## Output

Each script writes to `BLUEs_results/<LOC>/`:

- `trait_BLUEs_results.csv` — columns: `trait, location, accession,
  SampleName, BLUE, SE, Convergence_Status, Iterations, QQ_Correlation,
  Normality_Status`. One row per (trait × accession).
- `trait_heritability_results.csv` — per-trait generalised heritability +
  convergence flag.
- `diagnostic_plots/<trait>_diagnostics.pdf` — residuals vs fitted, Q-Q plot,
  residual histogram, response-by-year boxplot.

The two `trait_BLUEs_results.csv` files share the same schema and can be
concatenated with `rbind` for downstream phenotype plots.

## Usage

Run from `preprocessing/BLUEs/` (paths are relative):

```bash
# Phenotype-summary BLUEs (one value per accession per location)
Rscript BLUEs_AUS.R
Rscript BLUEs_PAK.R

# Downstream merge for phenotype-correlation plots
Rscript merge_BLUEs_for_phenotype_plots.R      # → phenotype-plot input

# GP model input (raw means everywhere — does not depend on BLUEs above)
cd .. && Rscript build_GP_phenotype_input.R
```

Convergence uses a retry loop (up to `max_iterations = 20`). Normality of
residuals is checked via Q-Q correlation (Shapiro-Wilk is over-sensitive at
this n); the interpretation string is written into the output CSV.
