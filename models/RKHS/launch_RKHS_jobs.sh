#!/bin/bash
################################################################################
# launch_RKHS_jobs.sh
#
# Generates and (optionally) submits SLURM job scripts for parallelised
# RKHS (multi-kernel) AUSPAK cross-validation.
#
# Job breakdown per trait:
#   - 1 job for CV1 (all 15 iterations × 5 folds)
#   - 1 job for CV2 (all 15 iterations)
#   - 1 job for CV0 + CrossLoc (deterministic, combined)
#   Total: 3 jobs per trait × 7 traits = 21 jobs
#
# RKHS is faster per fit than BayesC (kernel operations vs 400k+ marker
# effects), so all iterations run within a single job per CV scheme.
# Each script writes intermediate checkpoints, so interrupted runs are
# recoverable.
#
# Usage:
#   bash launch_RKHS_jobs.sh                    # generate scripts only
#   bash launch_RKHS_jobs.sh submit             # generate and submit all
#   bash launch_RKHS_jobs.sh submit DTF         # submit one trait only
#
# After all jobs finish, run the aggregation script:
#   Rscript aggregate_RKHS_results.R all
################################################################################

# ── SLURM parameters ────────────────────────────────────────────

PARTITION="batch"

# CV1: 15 iters × 5 folds = 75 BGLR fits (kernel-based, faster than BayesC)
TIME_CV1="110:00:00"
MEM_CV1="80G"                # kernels are ~551×551, much smaller than 400k markers

# CV2: 15 iterations × 5 folds = 75 BGLR fits (stratified CV, matching BayesC)
TIME_CV2="110:00:00"
MEM_CV2="80G"

# CV0 + CrossLoc: ~6 CV0 fits + 2 CrossLoc fits
TIME_DETERM="110:00:00"
MEM_DETERM="80G"

CPUS=1                       # BGLR Gibbs sampler is single-threaded

# MCMC settings (passed as environment variables to R scripts)
# Defaults match BayesC: 15000/5000/5. Override here if needed.
RKHS_NITER=""                # leave empty to use script default (15000)
RKHS_BURNIN=""               # leave empty to use script default (5000)
RKHS_THIN=""                 # leave empty to use script default (5)

ACCOUNT=""                   # leave empty if not needed
MAIL="stansccs@kaust.edu.sa"

# ── Paths  ───────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WORK_DIR="${SCRIPT_DIR}"     # where phenotype + genotype + R scripts live
MARKER_FILE="auspak_for_rkhs.raw"

# ── Traits ───────────────────────────────────────────────────────────────────

TRAITS=("DTF" "DTH" "PtHt" "PcleLng" "SdLen" "TGW" "SdW_z")

# ── Parse arguments ──────────────────────────────────────────────────────────

DO_SUBMIT=false
SINGLE_TRAIT=""

if [ "${1}" == "submit" ]; then
  DO_SUBMIT=true
  if [ -n "${2}" ]; then
    SINGLE_TRAIT="${2}"
    # Validate trait
    VALID=false
    for T in "${TRAITS[@]}"; do
      if [ "${T}" == "${SINGLE_TRAIT}" ]; then VALID=true; break; fi
    done
    if [ "${VALID}" == "false" ]; then
      echo "ERROR: Invalid trait '${SINGLE_TRAIT}'"
      echo "Valid traits: ${TRAITS[*]}"
      exit 1
    fi
    TRAITS=("${SINGLE_TRAIT}")
  fi
fi

# ── Helper: write common SBATCH header ───────────────────────────────────────

write_sbatch_header() {
  local JOB_SCRIPT=$1
  local JOB_NAME=$2
  local LOG_PREFIX=$3
  local TIME=$4
  local MEM=$5

  cat > "${JOB_SCRIPT}" <<EOF
#!/bin/bash
#SBATCH --job-name=${JOB_NAME}
#SBATCH --output=logs/${LOG_PREFIX}_%j.out
#SBATCH --error=logs/${LOG_PREFIX}_%j.err
#SBATCH --partition=${PARTITION}
#SBATCH --time=${TIME}
#SBATCH --mem=${MEM}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --nodes=1
#SBATCH --ntasks=1
EOF

  if [ -n "${ACCOUNT}" ]; then
    echo "#SBATCH --account=${ACCOUNT}" >> "${JOB_SCRIPT}"
  fi
  if [ -n "${MAIL}" ]; then
    echo "#SBATCH --mail-user=${MAIL}" >> "${JOB_SCRIPT}"
    echo "#SBATCH --mail-type=FAIL" >> "${JOB_SCRIPT}"
  fi
}

# ── Generate job scripts ─────────────────────────────────────────────────────

mkdir -p slurm_scripts logs
TOTAL_JOBS=0

for TRAIT in "${TRAITS[@]}"; do

  echo "=== Generating jobs for ${TRAIT} ==="

  # ── CV1: single job for all iterations ───────────────────────────────────
  JOB_SCRIPT="slurm_scripts/rkhs_cv1_${TRAIT}.slurm"

  write_sbatch_header "${JOB_SCRIPT}" \
    "RK1_${TRAIT}" \
    "rkhs_cv1_${TRAIT}" \
    "${TIME_CV1}" "${MEM_CV1}"

  cat >> "${JOB_SCRIPT}" <<EOF

echo "Job started: \$(date)"
echo "Node: \$(hostname) | Trait: ${TRAIT} | RKHS CV1 (all iterations)"

module load R/4.5.0
cd ${WORK_DIR}

# MCMC settings (empty = use script defaults)
${RKHS_NITER:+export RKHS_NITER=${RKHS_NITER}}
${RKHS_BURNIN:+export RKHS_BURNIN=${RKHS_BURNIN}}
${RKHS_THIN:+export RKHS_THIN=${RKHS_THIN}}

Rscript RKHS_CV1.R ${TRAIT} ${MARKER_FILE}

echo "Job finished: \$(date)"
EOF
  chmod +x "${JOB_SCRIPT}"
  TOTAL_JOBS=$((TOTAL_JOBS + 1))
  echo "  CV1: 1 job script"

  # ── CV2: single job for all iterations ───────────────────────────────────
  JOB_SCRIPT="slurm_scripts/rkhs_cv2_${TRAIT}.slurm"

  write_sbatch_header "${JOB_SCRIPT}" \
    "RK2_${TRAIT}" \
    "rkhs_cv2_${TRAIT}" \
    "${TIME_CV2}" "${MEM_CV2}"

  cat >> "${JOB_SCRIPT}" <<EOF

echo "Job started: \$(date)"
echo "Node: \$(hostname) | Trait: ${TRAIT} | RKHS CV2 (all iterations)"

module load R/4.5.0
cd ${WORK_DIR}

# MCMC settings (empty = use script defaults)
${RKHS_NITER:+export RKHS_NITER=${RKHS_NITER}}
${RKHS_BURNIN:+export RKHS_BURNIN=${RKHS_BURNIN}}
${RKHS_THIN:+export RKHS_THIN=${RKHS_THIN}}

Rscript RKHS_CV2.R ${TRAIT} ${MARKER_FILE}

echo "Job finished: \$(date)"
EOF
  chmod +x "${JOB_SCRIPT}"
  TOTAL_JOBS=$((TOTAL_JOBS + 1))
  echo "  CV2: 1 job script"

  # ── CV0 + CrossLoc: single combined job ──────────────────────────────────
  JOB_SCRIPT="slurm_scripts/rkhs_cv0xl_${TRAIT}.slurm"

  write_sbatch_header "${JOB_SCRIPT}" \
    "RK0X_${TRAIT}" \
    "rkhs_cv0xl_${TRAIT}" \
    "${TIME_DETERM}" "${MEM_DETERM}"

  cat >> "${JOB_SCRIPT}" <<EOF

echo "Job started: \$(date)"
echo "Node: \$(hostname) | Trait: ${TRAIT} | RKHS CV0 + CrossLoc"

module load R/4.5.0
cd ${WORK_DIR}

# MCMC settings (empty = use script defaults)
${RKHS_NITER:+export RKHS_NITER=${RKHS_NITER}}
${RKHS_BURNIN:+export RKHS_BURNIN=${RKHS_BURNIN}}
${RKHS_THIN:+export RKHS_THIN=${RKHS_THIN}}

Rscript RKHS_CV0_CrossLoc.R ${TRAIT} ${MARKER_FILE}

echo "Job finished: \$(date)"
EOF
  chmod +x "${JOB_SCRIPT}"
  TOTAL_JOBS=$((TOTAL_JOBS + 1))
  echo "  CV0+CrossLoc: 1 job script"

  echo "  Subtotal for ${TRAIT}: 3 jobs"
  echo ""
done

echo "========================================"
echo "Generated ${TOTAL_JOBS} job scripts in slurm_scripts/"
echo "========================================"

# ── Submit if requested ──────────────────────────────────────────────────────

if [ "${DO_SUBMIT}" == "true" ]; then
  echo ""
  echo "Submitting jobs..."
  SUBMITTED=0

  for TRAIT in "${TRAITS[@]}"; do
    echo "--- ${TRAIT} ---"

    # Submit CV1
    JOB_ID=$(sbatch "slurm_scripts/rkhs_cv1_${TRAIT}.slurm" | awk '{print $NF}')
    echo "  CV1: job ${JOB_ID}"
    SUBMITTED=$((SUBMITTED + 1))

    # Submit CV2
    JOB_ID=$(sbatch "slurm_scripts/rkhs_cv2_${TRAIT}.slurm" | awk '{print $NF}')
    echo "  CV2: job ${JOB_ID}"
    SUBMITTED=$((SUBMITTED + 1))

    # Submit CV0+CrossLoc
    JOB_ID=$(sbatch "slurm_scripts/rkhs_cv0xl_${TRAIT}.slurm" | awk '{print $NF}')
    echo "  CV0+CrossLoc: job ${JOB_ID}"
    SUBMITTED=$((SUBMITTED + 1))
  done

  echo ""
  echo "Submitted ${SUBMITTED} jobs. Monitor with: squeue -u \$USER"
  echo ""
  echo "After all jobs finish, aggregate results with:"
  echo "  Rscript aggregate_RKHS_results.R all"
else
  echo ""
  echo "To submit all traits:     bash launch_RKHS_jobs.sh submit"
  echo "To submit one trait:      bash launch_RKHS_jobs.sh submit DTF_blue"
  echo ""
  echo "After all jobs finish:    Rscript aggregate_RKHS_results.R all"
fi
