#!/bin/bash

# Submit calibrated global-unweighted AAO radiative mode-3 production for
# the padded RGA 10.604-GeV analysis phase space.
#
# Target:
#   100M generated events = 20,000 jobs * 5,000 events/job
#
# Mode 3 samples a 75/25 mixture of
#   direct: (1/Q_l^2, x_l, radiative variables, -t_h, phi_h)
#   legacy: (1/Q_l^2, E',  radiative variables, cos(theta*), phi_h)
# and emits globally equal-weight LUND events in the padded final-state box.
# A compatible envelope calibration is mandatory.
#
# Usage:
#   ./submit_aao_rad_mode3_rga10604.sh /path/to/envelope_calibration.json
#   ./submit_aao_rad_mode3_rga10604.sh --prepare-only /path/to/envelope_calibration.json

set -euo pipefail

PREPARE_ONLY=0
if [ "${1:-}" = "--prepare-only" ]; then
    PREPARE_ONLY=1
    shift
fi
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 [--prepare-only] /path/to/envelope_calibration.json" >&2
    exit 2
fi

# ---------------- CONFIG ----------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
AAO_RAD_DIR="${SCRIPT_DIR}/aao_rad"
DRIVER="${AAO_RAD_DIR}/radiative_mode3.py"
EXECUTABLE="${AAO_RAD_DIR}/build/aao_rad"
INPUT_TEMPLATE="${AAO_RAD_DIR}/input/rga10604_mode3_legacy.inp"
CONFIG="${REPO_DIR}/configs/analysis/rga/10.604.json"
CALIBRATION="$(realpath "$1")"

OUTBASE="/volatile/clas12/storyf/rad_mode3"
SCRATCH_ROOT="/scratch/${USER}/tmpfs"

BEAM="10.604"
Q2MIN="1.0"
Q2MAX="10.5"
EPMIN="2.0"
EPMAX="10.604"
EGMIN="0.010"
PADDING="0.035"
DIRECT_FRACTION="0.75"

LABEL="stage2_mode3"
NEVENTS_PER_JOB=5000
NJOBS=20000
TOTAL_EVENTS=$((NEVENTS_PER_JOB * NJOBS))
FILES_PER_DIRECTORY=5000
JOBS_PER_WORKFLOW=5000
SEED_BASE=1307000001
WALLTIME="24hr"
RAM="2gb"
DISK="2gb"
# ----------------------------------------

for required in "$DRIVER" "$EXECUTABLE" "$INPUT_TEMPLATE" "$CONFIG" "$CALIBRATION"; do
    if [ ! -f "$required" ]; then
        echo "ERROR: required file not found: $required" >&2
        exit 1
    fi
done
if [ ! -x "$EXECUTABLE" ]; then
    echo "ERROR: mode-3 executable is not executable: $EXECUTABLE" >&2
    exit 1
fi
if [ $((NJOBS % JOBS_PER_WORKFLOW)) -ne 0 ]; then
    echo "ERROR: NJOBS must be divisible by JOBS_PER_WORKFLOW" >&2
    exit 1
fi

BTAG="E${BEAM}"
Q2TAG="Q2_${Q2MIN}_${Q2MAX}"
EPTAG="EP_${EPMIN}_${EPMAX}"
EGTAG="EG_${EGMIN}"
PADTAG="P_${PADDING}"
WORKFLOW_BASE="aao_rad_${LABEL}_${BTAG}_${Q2TAG}_${EPTAG}_${EGTAG}_${PADTAG}"
DIR="${OUTBASE}/${WORKFLOW_BASE}"
REVISION="$(git -C "$SCRIPT_DIR" rev-parse HEAD)"

echo "Preparing mode-3 production: ${WORKFLOW_BASE}"
echo "  nominal analysis: Q2 ${Q2MIN}-${Q2MAX}, EP >= ${EPMIN} GeV"
echo "  padding fraction: ${PADDING} in Q2, xB, and -t"
echo "  minimum photon energy: ${EGMIN} GeV"
echo "  direct/legacy proposal fractions: ${DIRECT_FRACTION}/0.25"
echo "  events/job=${NEVENTS_PER_JOB}"
echo "  jobs=${NJOBS} (four SWIF workflows of ${JOBS_PER_WORKFLOW})"
echo "  total events=${TOTAL_EVENTS}"
echo "  output=${DIR}"
echo "  calibration=${CALIBRATION}"
echo "  generator revision=${REVISION}"

python3 "$DRIVER" prepare-production \
    --config "$CONFIG" \
    --input "$INPUT_TEMPLATE" \
    --calibration "$CALIBRATION" \
    --output "$DIR" \
    --tag "$WORKFLOW_BASE" \
    --padding-fraction "$PADDING" \
    --direct-fraction "$DIRECT_FRACTION" \
    --minimum-photon-energy "$EGMIN" \
    --total-events "$TOTAL_EVENTS" \
    --events-per-job "$NEVENTS_PER_JOB" \
    --lund-files-per-directory "$FILES_PER_DIRECTORY" \
    --heartbeat-interval 100000 \
    --seed-base "$SEED_BASE" \
    --generator-revision "$REVISION"

submit_scripts=()
for start in $(seq 0 "$JOBS_PER_WORKFLOW" $((NJOBS - JOBS_PER_WORKFLOW))); do
    stop=$((start + JOBS_PER_WORKFLOW))
    chunk=$((start / JOBS_PER_WORKFLOW))
    workflow="${WORKFLOW_BASE}_c$(printf '%02d' "$chunk")"
    submit_script="${DIR}/submit_${workflow}.sh"

    python3 "$DRIVER" emit-swif \
        "${DIR}/manifest.json" \
        --workflow "$workflow" \
        --executable "$EXECUTABLE" \
        --task-start "$start" \
        --task-stop "$stop" \
        --cores 1 \
        --ram "$RAM" \
        --disk "$DISK" \
        --walltime "$WALLTIME" \
        --scratch-root "$SCRATCH_ROOT" \
        --output "$submit_script"

    submit_scripts+=("$submit_script")
done

echo
echo "Prepared ${#submit_scripts[@]} SWIF submission scripts:"
printf '  %s\n' "${submit_scripts[@]}"

if [ "$PREPARE_ONLY" -eq 1 ]; then
    echo
    echo "Prepare-only requested; no workflows were submitted."
    echo "Inspect the manifest, then run each submit script listed above."
    exit 0
fi

for submit_script in "${submit_scripts[@]}"; do
    "$submit_script"
done

echo
echo "Submitted all four mode-3 workflow chunks."
echo "Check with:"
for chunk in 0 1 2 3; do
    workflow="${WORKFLOW_BASE}_c$(printf '%02d' "$chunk")"
    echo "  swif2 status ${workflow} -summary -problems"
done
echo
echo "After all jobs finish:"
echo "  python3 ${DRIVER} status ${DIR}/manifest.json"
echo "  python3 ${DRIVER} finalize-production ${DIR}/manifest.json"
