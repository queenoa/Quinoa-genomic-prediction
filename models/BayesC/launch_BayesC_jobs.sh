#!/bin/bash
################################################################################
# launch_BayesC_jobs.sh
#
# Generates and (optionally) submits SLURM job scripts for parallelised
# BayesC AUSPAK cross-validation.
#
# Job breakdown per trait (full run):
#   - 15 jobs for CV1 (one per iteration, each loops over 5 folds)
#   - 15 jobs for CV2 (one per iteration, each loops over 5 folds)
#   - 1  job for CV0 (leave-one-location-year-out)
#   - 1  job for CrossLoc (cross-location transfer)
#   Total: 33 jobs per trait × 7 traits = 231 jobs
#
# Usage:
#   bash launch_BayesC_jobs.sh                           # generate scripts only (all iters)
#   bash launch_BayesC_jobs.sh submit                    # generate and submit all
#   bash launch_BayesC_jobs.sh submit 1 2                # all traits, iters 1-2 only
#   bash launch_BayesC_jobs.sh submit cv0                # all traits, CV0 only
#   bash launch_BayesC_jobs.sh submit crossloc           # all traits, CrossLoc only
#   bash launch_BayesC_jobs.sh submit DTF_blue           # one trait, all jobs
#   bash launch_BayesC_jobs.sh submit DTF_blue 1 2       # one trait, iters 1-2 only
#   bash launch_BayesC_jobs.sh submit DTF_blue 3 5       # one trait, resume iters 3-5
#   bash launch_BayesC_jobs.sh submit DTF_blue cv0       # one trait, CV0 only
#   bash launch_BayesC_jobs.sh submit DTF_blue crossloc  # one trait, CrossLoc only
#   bash launch_BayesC_jobs.sh submit fullmodel           # all traits, full model only
#   bash launch_BayesC_jobs.sh submit DTF_blue fullmodel  # one trait, full model only
#
# After all jobs finish, run the aggregation script:
#   Rscript aggregate_BayesC_results.R all
################################################################################

# ── SLURM parameters (EDIT THESE) ────────────────────────────────────────────

# CV1 and CV2 per-iteration jobs (both loop over 5 folds)
PARTITION="batch"
TIME_ITER="200:00:00"         # per-iteration wall time (5 BGLR fits)
MEM_ITER="80G"               # ~400k markers × ~1900 obs
CPUS=1                       # BGLR Gibbs sampler is single-threaded

# CV0, CrossLoc, and full model (separate jobs, each deterministic single-pass)
TIME_CV0="200:00:00"         # ~6 leave-one-location-year-out fits
TIME_CROSSLOC="200:00:00"   # ~2 cross-location fits (fewer obs per fit)
TIME_FULLMODEL="200:00:00"  # single fit on all data
MEM_DETERM="80G"

ACCOUNT=""                   # leave empty if not needed
MAIL=""  

# ── Paths (EDIT THESE) ───────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WORK_DIR="${SCRIPT_DIR}"     # where phenotype + genotype + R scripts live
MARKER_FILE="pruned05_AUSPAK_for_bayesC.raw"

# ── Traits and iterations ────────────────────────────────────────────────────

TRAITS=("DTF_blue" "DTH_blue" "PtHt_blue" "PcleLng_blue" "SdLen_blue" "TGW_blue" "SdW_z_blue")
N_ITERATIONS=15

# ── Parse arguments ──────────────────────────────────────────────────────────

DO_SUBMIT=false
START_ITER=1
END_ITER=${N_ITERATIONS}
RUN_CV1_CV2=true
RUN_CV0=true
RUN_CROSSLOC=true
RUN_FULLMODEL=false

# Helper: parse mode keyword or iteration range from positional args
parse_mode_args() {
  local ARG1=$1
  local ARG2=$2
  if [ -n "${ARG1}" ]; then
    if [ "${ARG1}" == "cv0" ]; then
      RUN_CV1_CV2=false
      RUN_CROSSLOC=false
    elif [ "${ARG1}" == "crossloc" ]; then
      RUN_CV1_CV2=false
      RUN_CV0=false
    elif [ "${ARG1}" == "fullmodel" ]; then
      RUN_CV1_CV2=false
      RUN_CV0=false
      RUN_CROSSLOC=false
      RUN_FULLMODEL=true
    else
      START_ITER=${ARG1}
      END_ITER=${ARG2:-${ARG1}}
      RUN_CV0=false
      RUN_CROSSLOC=false
      if [ "${START_ITER}" -lt 1 ] || [ "${END_ITER}" -gt "${N_ITERATIONS}" ] || [ "${START_ITER}" -gt "${END_ITER}" ]; then
        echo "ERROR: Invalid iteration range ${START_ITER}-${END_ITER} (N_ITERATIONS=${N_ITERATIONS})"
        exit 1
      fi
    fi
  fi
}

if [ "${1}" == "submit" ]; then
  DO_SUBMIT=true
  if [ -n "${2}" ]; then
    # Check if arg 2 is a trait name, a number, or a mode keyword
    IS_TRAIT=false
    for T in "${TRAITS[@]}"; do
      if [ "${T}" == "${2}" ]; then IS_TRAIT=true; break; fi
    done

    if [ "${IS_TRAIT}" == "true" ]; then
      # submit <trait> [iter_start [iter_end] | cv0 | crossloc | fullmodel]
      TRAITS=("${2}")
      parse_mode_args "${3}" "${4}"
    elif [ "${2}" == "cv0" ] || [ "${2}" == "crossloc" ] || [ "${2}" == "fullmodel" ] || [[ "${2}" =~ ^[0-9]+$ ]]; then
      # submit <mode> | submit <iter_start> [iter_end]  (all traits)
      parse_mode_args "${2}" "${3}"
    else
      echo "ERROR: Invalid argument '${2}'"
      echo "Expected a trait name, iteration number, 'cv0', 'crossloc', or 'fullmodel'"
      echo "Valid traits: ${TRAITS[*]}"
      exit 1
    fi
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
    echo "#SBATCH --mail-type=END,FAIL" >> "${JOB_SCRIPT}"
  fi
}

# ── Generate job scripts ─────────────────────────────────────────────────────

mkdir -p slurm_scripts logs
TOTAL_JOBS=0

for TRAIT in "${TRAITS[@]}"; do
  
  echo "=== Generating jobs for ${TRAIT} ==="

  # ── CV1: one job per iteration ───────────────────────────────────────────
  if [ "${RUN_CV1_CV2}" == "true" ]; then
  for ITER in $(seq ${START_ITER} ${END_ITER}); do
    ITER_PAD=$(printf "%02d" ${ITER})
    JOB_SCRIPT="slurm_scripts/bc_cv1_${TRAIT}_i${ITER_PAD}.slurm"

    write_sbatch_header "${JOB_SCRIPT}" \
      "BC1_${TRAIT}_i${ITER_PAD}" \
      "bc_cv1_${TRAIT}_i${ITER_PAD}" \
      "${TIME_ITER}" "${MEM_ITER}"

    cat >> "${JOB_SCRIPT}" <<EOF

echo "Job started: \$(date)"
echo "Node: \$(hostname) | Trait: ${TRAIT} | CV1 iteration: ${ITER}"

module load R/4.5.0
cd ${WORK_DIR}

Rscript BayesC_CV1_single_iter.R ${TRAIT} ${ITER} ${MARKER_FILE}

echo "Job finished: \$(date)"
EOF
    chmod +x "${JOB_SCRIPT}"
    TOTAL_JOBS=$((TOTAL_JOBS + 1))
  done
  echo "  CV1: $((END_ITER - START_ITER + 1)) job scripts (iters ${START_ITER}-${END_ITER})"

  # ── CV2: one job per iteration ───────────────────────────────────────────
  for ITER in $(seq ${START_ITER} ${END_ITER}); do
    ITER_PAD=$(printf "%02d" ${ITER})
    JOB_SCRIPT="slurm_scripts/bc_cv2_${TRAIT}_i${ITER_PAD}.slurm"

    write_sbatch_header "${JOB_SCRIPT}" \
      "BC2_${TRAIT}_i${ITER_PAD}" \
      "bc_cv2_${TRAIT}_i${ITER_PAD}" \
      "${TIME_ITER}" "${MEM_ITER}"

    cat >> "${JOB_SCRIPT}" <<EOF

echo "Job started: \$(date)"
echo "Node: \$(hostname) | Trait: ${TRAIT} | CV2 iteration: ${ITER}"

module load R/4.5.0
cd ${WORK_DIR}

Rscript BayesC_CV2_single_iter.R ${TRAIT} ${ITER} ${MARKER_FILE}

echo "Job finished: \$(date)"
EOF
    chmod +x "${JOB_SCRIPT}"
    TOTAL_JOBS=$((TOTAL_JOBS + 1))
  done
  echo "  CV2: $((END_ITER - START_ITER + 1)) job scripts (iters ${START_ITER}-${END_ITER})"
  fi

  # ── CV0: leave-one-location-year-out ──────────────────────────────────────
  if [ "${RUN_CV0}" == "true" ]; then
  JOB_SCRIPT="slurm_scripts/bc_cv0_${TRAIT}.slurm"

  write_sbatch_header "${JOB_SCRIPT}" \
    "BC0_${TRAIT}" \
    "bc_cv0_${TRAIT}" \
    "${TIME_CV0}" "${MEM_DETERM}"

  cat >> "${JOB_SCRIPT}" <<EOF

echo "Job started: \$(date)"
echo "Node: \$(hostname) | Trait: ${TRAIT} | CV0 (leave-one-location-year-out)"

module load R/4.5.0
cd ${WORK_DIR}

Rscript BayesC_CV0.R ${TRAIT} ${MARKER_FILE}

echo "Job finished: \$(date)"
EOF
  chmod +x "${JOB_SCRIPT}"
  TOTAL_JOBS=$((TOTAL_JOBS + 1))
  echo "  CV0: 1 job script"
  fi

  # ── CrossLoc: cross-location transfer ──────────────────────────────────────
  if [ "${RUN_CROSSLOC}" == "true" ]; then
  JOB_SCRIPT="slurm_scripts/bc_crossloc_${TRAIT}.slurm"

  write_sbatch_header "${JOB_SCRIPT}" \
    "BCX_${TRAIT}" \
    "bc_crossloc_${TRAIT}" \
    "${TIME_CROSSLOC}" "${MEM_DETERM}"

  cat >> "${JOB_SCRIPT}" <<EOF

echo "Job started: \$(date)"
echo "Node: \$(hostname) | Trait: ${TRAIT} | Cross-location transfer"

module load R/4.5.0
cd ${WORK_DIR}

Rscript BayesC_CrossLoc.R ${TRAIT} ${MARKER_FILE}

echo "Job finished: \$(date)"
EOF
  chmod +x "${JOB_SCRIPT}"
  TOTAL_JOBS=$((TOTAL_JOBS + 1))
  echo "  CrossLoc: 1 job script"
  fi

  # ── Full model: train on all data ──────────────────────────────────────────
  if [ "${RUN_FULLMODEL}" == "true" ]; then
  JOB_SCRIPT="slurm_scripts/bc_fullmodel_${TRAIT}.slurm"

  write_sbatch_header "${JOB_SCRIPT}" \
    "BCF_${TRAIT}" \
    "bc_fullmodel_${TRAIT}" \
    "${TIME_FULLMODEL}" "${MEM_DETERM}"

  cat >> "${JOB_SCRIPT}" <<EOF

echo "Job started: \$(date)"
echo "Node: \$(hostname) | Trait: ${TRAIT} | Full model (all data)"

module load R/4.5.0
cd ${WORK_DIR}

Rscript BayesC_full_model.R ${TRAIT} ${MARKER_FILE}

echo "Job finished: \$(date)"
EOF
  chmod +x "${JOB_SCRIPT}"
  TOTAL_JOBS=$((TOTAL_JOBS + 1))
  echo "  Full model: 1 job script"
  fi

  echo "  Subtotal for ${TRAIT}: ${TOTAL_JOBS} jobs"
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

    if [ "${RUN_CV1_CV2}" == "true" ]; then
    # Submit CV1 iterations
    for ITER in $(seq ${START_ITER} ${END_ITER}); do
      ITER_PAD=$(printf "%02d" ${ITER})
      JOB_ID=$(sbatch "slurm_scripts/bc_cv1_${TRAIT}_i${ITER_PAD}.slurm" | awk '{print $NF}')
      echo "  CV1 iter ${ITER}: job ${JOB_ID}"
      SUBMITTED=$((SUBMITTED + 1))
    done

    # Submit CV2 iterations
    for ITER in $(seq ${START_ITER} ${END_ITER}); do
      ITER_PAD=$(printf "%02d" ${ITER})
      JOB_ID=$(sbatch "slurm_scripts/bc_cv2_${TRAIT}_i${ITER_PAD}.slurm" | awk '{print $NF}')
      echo "  CV2 iter ${ITER}: job ${JOB_ID}"
      SUBMITTED=$((SUBMITTED + 1))
    done
    fi

    if [ "${RUN_CV0}" == "true" ]; then
    # Submit CV0
    JOB_ID=$(sbatch "slurm_scripts/bc_cv0_${TRAIT}.slurm" | awk '{print $NF}')
    echo "  CV0: job ${JOB_ID}"
    SUBMITTED=$((SUBMITTED + 1))
    fi

    if [ "${RUN_CROSSLOC}" == "true" ]; then
    # Submit CrossLoc
    JOB_ID=$(sbatch "slurm_scripts/bc_crossloc_${TRAIT}.slurm" | awk '{print $NF}')
    echo "  CrossLoc: job ${JOB_ID}"
    SUBMITTED=$((SUBMITTED + 1))
    fi

    if [ "${RUN_FULLMODEL}" == "true" ]; then
    # Submit full model
    JOB_ID=$(sbatch "slurm_scripts/bc_fullmodel_${TRAIT}.slurm" | awk '{print $NF}')
    echo "  Full model: job ${JOB_ID}"
    SUBMITTED=$((SUBMITTED + 1))
    fi
  done
  
  echo ""
  echo "Submitted ${SUBMITTED} jobs. Monitor with: squeue -u \$USER"
  echo ""
  echo "After all jobs finish, aggregate results with:"
  echo "  Rscript aggregate_BayesC_results.R all"
else
  echo ""
  echo "To submit all traits:     bash launch_BayesC_jobs.sh submit"
  echo "To submit one trait:      bash launch_BayesC_jobs.sh submit DTF_blue"
  echo ""
  echo "After all jobs finish:    Rscript aggregate_BayesC_results.R all"
fi