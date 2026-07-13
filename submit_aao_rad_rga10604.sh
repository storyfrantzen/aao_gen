#!/bin/bash

# Submit AAO radiative pi0 production for the RGA 10.604 GeV analysis.
#
# RGA analysis binning:
#   Q2: 1.0-10.5 GeV^2
#   xB: 0.05-0.70
#   -t: 0.09-2.0 GeV^2
#   phi: 0-360 deg in 18 deg bins
#
# The AAO radiative generator input controls Q2 and scattered-electron energy
# directly.  The default EP lower bound follows the analysis y_max=0.8 cut:
#   E' = Ebeam * (1 - y_max) = 10.604 * 0.2 = 2.1208 GeV.

set -euo pipefail

# ---------------- CONFIG ----------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="${SCRIPT_DIR}/run_swif_job.sh"
OUTBASE="${OUTBASE:-/volatile/clas12/${USER}/LUND/aao_rad_production}"

BEAM="${BEAM:-10.604}"
Q2MIN="${Q2MIN:-1.0}"
Q2MAX="${Q2MAX:-10.5}"
EPMIN="${EPMIN:-2.10}"
EPMAX="${EPMAX:-10.0}"
EGMIN="${EGMIN:-0.010}"

LABEL="${LABEL:-rga10604_stage2_rad}"
NEVENTS_PER_JOB="${NEVENTS_PER_JOB:-5000}"
NJOBS="${NJOBS:-20000}"
WALLTIME="${WALLTIME:-24hr}"
RAM="${RAM:-256mb}"
DISK="${DISK:-1gb}"
# ----------------------------------------

if [ ! -x "$SCRIPT" ]; then
    echo "ERROR: SWIF job script is not executable or not found: $SCRIPT" >&2
    exit 1
fi

BTAG="E${BEAM}"
Q2TAG="Q2_${Q2MIN}_${Q2MAX}"
EPTAG="EP_${EPMIN}_${EPMAX}"
EGTAG="EG_${EGMIN}"

WORKFLOW="aao_rad_${LABEL}_${BTAG}_${Q2TAG}_${EPTAG}_${EGTAG}"
DIR="${OUTBASE}/${WORKFLOW}"
INPUT_DIR="${DIR}/inputs"

mkdir -p "$INPUT_DIR"

echo "Creating workflow: ${WORKFLOW}"
echo "  phase space: Q2 ${Q2MIN}-${Q2MAX}, EP ${EPMIN}-${EPMAX}, EG ${EGMIN}"
echo "  events/job=${NEVENTS_PER_JOB}"
echo "  jobs=${NJOBS}"
echo "  total events=$((NEVENTS_PER_JOB * NJOBS))"
echo "  walltime=${WALLTIME}"
echo "  ram=${RAM}"
echo "  disk=${DISK}"
echo "  output=${DIR}"

swif2 create -workflow "$WORKFLOW"

for job in $(seq 1 "$NJOBS"); do
    JOB_ID=$(printf "%05d" "$job")
    INPUT="${INPUT_DIR}/input_${JOB_ID}.inp"

    # Input structure expected by aaorad_gen:
    #  1 th_opt
    #  2 flag_ehel
    #  3 reg1 reg2 reg3 reg4
    #  4 npart
    #  5 epirea
    #  6 mm_cut
    #  7 target thickness
    #  8 target radius
    #  9 vertex_x
    # 10 vertex_y
    # 11 vz
    # 12 ebeam
    # 13 q2_min q2_max
    # 14 ep_min ep_max
    # 15 delta / minimum photon energy
    # 16 nmax
    # 17 fmcall
    cat > "$INPUT" << EOF
5
1
.20 .12 .20 .20
4
1
.2
5.0
0.8
0.0
0.0
0.0
${BEAM}
${Q2MIN} ${Q2MAX}
${EPMIN} ${EPMAX}
${EGMIN}
${NEVENTS_PER_JOB}
1
EOF

    swif2 add-job \
        -workflow "$WORKFLOW" \
        -name "${WORKFLOW}_${JOB_ID}" \
        -cores 1 \
        -disk "$DISK" \
        -ram "$RAM" \
        -time "$WALLTIME" \
        -os el9 \
        -input "$INPUT" "$INPUT" \
        -- "$SCRIPT" "$INPUT" "$DIR"
done

swif2 run "$WORKFLOW"

echo "Submitted ${WORKFLOW}"
echo
echo "Check with:"
echo "  swif2 diagnose ${WORKFLOW}"
echo "  swif2 status ${WORKFLOW}"
echo "  swif2 problems ${WORKFLOW}"
