# RKHS Pipeline

BGLR RKHS (Reproducing Kernel Hilbert Space) kernel-averaging implementation for the AUSPAK quinoa multi-environment genomic prediction comparison. See the repository [README](../../README.md) for the project overview, trait/environment definitions, CV-scheme designs, and evaluation metrics shared across all four pipelines.

Where BayesC models explicit SNP marker effects, RKHS captures genomic relationships through Gaussian kernels computed from the full marker matrix. Both pipelines share the same fixed-effect structure, preprocessing, evaluation metrics, and CV designs.

## RKHS method

### Gaussian kernel construction

1. Load PLINK `.raw` marker matrix, mean-impute missing genotypes
2. Center and scale each SNP (zero mean, unit variance) so all markers contribute equally regardless of MAF
3. Compute squared Euclidean distance matrix: `D[i,j] = ||x_i - x_j||^2`
4. Normalise: `D <- D / mean(D)`
5. Compute three Gaussian kernels at different bandwidths `K(x_i, x_j) = exp(-h * D[i,j])`:
   - `h1 = 1/(5 * median(D))` — wide (smooth, captures broad similarity)
   - `h2 = 1/median(D)` — medium
   - `h3 = 5/median(D)` — narrow (captures fine-grained similarity)

Following Perez & de los Campos (2014, Box 11) and Cuevas et al. (2016). BGLR estimates the variance component for each kernel separately, effectively performing Bayesian model averaging over bandwidths (kernel averaging, KA).

### Kernel expansion to observation level

For multi-environment data in long format, observation `i` corresponds to genotype `g_i`. The observation-level kernel is `K_obs[i,j] = K_geno[g_i, g_j]` — the kernel analogue of `build_obs_marker_matrix()` in BayesC, but much cheaper (indexing vs. matrix multiplication).

Kernels are computed once from the full marker matrix and checkpointed to `RKHS_kernels.RData`.

### Empty location-years are retained

Unlike BayesC, empty location-years (100% missing for the target trait) are **not** dropped. The kernel is precomputed at the genotype level and expanded to observation level by simple indexing, so the cost of extra NA rows is negligible.

## Model structure

Defined in `RKHS_utils.R`. Each BGLR fit includes:

- **Fixed effects**: location + year-within-location, via `build_fixed_design()` — identical to BayesC; all fixed since BGLR lacks random effects
- **3 RKHS kernels**: Gaussian kernels at narrow, medium, and wide bandwidths — BGLR estimates a variance component for each (Bayesian kernel averaging)
- **Residual**: single residual variance — BGLR does not support heterogeneous residuals (`groups`) with the RKHS model type, so unlike BayesC, RKHS uses one residual variance across all environments
- **Preprocessing**: z-score standardisation per location-year; mean-imputation of missing marker dosages

MCMC settings: 15,000 iterations, 5,000 burn-in, thinning every 5. Overridable via env vars `RKHS_NITER`, `RKHS_BURNIN`, `RKHS_THIN`.

### CrossLoc variant

For cross-location prediction:
- **No location fixed effect** (target location cannot contribute to estimating it)
- **Year fixed effect within training location** (if >1 year), with target rows set to the reference level
- **3 RKHS kernels** built from the combined (train + target) genotype set
- Target genotypes have `y = NA`, so BGLR predicts them from kernel similarity to training genotypes

This is the kernel analogue of BayesC CrossLoc (which extracts marker effects and computes GEBVs). Prediction flows through the kernel rather than explicit marker effects.

## Files

| File | Purpose |
|------|---------|
| `RKHS_utils.R` | Shared utilities: data loading, kernel computation/checkpointing, z-score scaling, observation-level kernel expansion, fixed-effect design matrix, BGLR RKHS fitting wrapper, evaluation metrics |
| `RKHS_CV1.R` | All 15 iterations of CV1 (new genotypes). One SLURM job per trait. |
| `RKHS_CV2.R` | All 15 iterations of CV2 (sparse testing, 5-fold stratified). One SLURM job per trait. |
| `RKHS_CV0_CrossLoc.R` | Both CV0 and CrossLoc for one trait (deterministic, incremental CSV saving) |
| `launch_RKHS_jobs.sh` | Generates and optionally submits SLURM jobs (3 per trait) |
| `aggregate_RKHS_results.R` | Collects per-scheme CSVs into combined results and summary statistics |

## Usage

```bash
# Generate SLURM scripts only
bash launch_RKHS_jobs.sh

# Generate and submit all traits
bash launch_RKHS_jobs.sh submit

# Submit one trait only
bash launch_RKHS_jobs.sh submit DTF

# After all jobs finish, aggregate
Rscript aggregate_RKHS_results.R all
Rscript aggregate_RKHS_results.R DTF
```

## Input files

- **Phenotype file**: `../../data/AUSPAK_phenotypes_GP_input.csv` — columns `sample.id`, `location`, `year`, `location_year`, plus one column per trait (per-location means)
- **Marker file**: `../../data/auspak_for_rkhs.raw` — PLINK `.raw` format with IID column and SNP dosages (first 6 columns are PLINK metadata)
- **Kernel checkpoint**: `RKHS_kernels.RData` — computed once from the marker file and reused across all traits and CV schemes

## Output files

Per trait, per scheme:
- `cv_results_{SCHEME}_{TRAIT}_RKHS.csv` — accuracy metrics per location-year
- `predictions_{SCHEME}_{TRAIT}_RKHS.csv` — individual-level predictions

After aggregation:
- `cv_results_{TRAIT}_RKHS_all_schemes.csv` — combined row-level results
- `cv_summary_{TRAIT}_RKHS_all_schemes.csv` — mean/SD/min/max per scheme × location
- `predictions_{TRAIT}_RKHS_all_schemes.csv` — combined individual predictions

## Functions shared with BayesC/GBLUP

Identical-logic helpers in `RKHS_utils.R` (mirrored in `BayesC_utils.R` and `GBLUP_utils.R`):

- `VALID_TRAITS`, `LOWER_IS_BETTER_TRAITS`
- `apply_location_year_scaling()` — z-score by location-year
- `build_fixed_design()` — location + year-within-location design matrix
- `build_groups()` — environment group vector (defined but not used in RKHS fitting; BGLR does not support `groups` with the RKHS model type)
- `calculate_ndcg()` — NDCG@k
- `evaluate_predictions()` — Pearson, Spearman, NDCG@10
- `evaluate_per_location_year()` — index-based (BGLR row-aligned); GBLUP uses a join-based variant

## References

- Perez P, de los Campos G (2014) Genome-wide regression and prediction with the BGLR statistical package. *Genetics* 198:483-495
- de los Campos G et al. (2010) Semi-parametric genomic-enabled prediction of genetic values. *Genetics Research* 92:295-308
- Cuevas J et al. (2016) Genomic prediction of genotype × environment interaction kernel regression models. *Plant Genome* 9(3)
