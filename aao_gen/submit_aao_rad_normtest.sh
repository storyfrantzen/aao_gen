#!/bin/bash

# Submit stage-2 AAO radiative production for the RGK 6.535 GeV analysis.
#
# Stage 2 target:
#   100M generated events = 5000 jobs * 20000 events/job
#
# Nominal production point:
#   Q2 0.7-5.0, EP 1.00-4.5, EG 0.010
#
# Timing basis:
#   The broad EG=0.010 diagnostic had a worst observed walltime of about
#   9h55 per 10k events, so 20k events/job with a 24h request keeps margin.

set -euo pipefail

# ---------------- CONFIG ----------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="${SCRIPT_DIR}/run_swif_job.s
OUTBASE="/volatile/clas12/storyf/LUND/aao_rad_normtest"

BEAM="6.535"
Q2MIN="0.7"
Q2MAX="5.0"
EPMIN="1.00"
EPMAX="4.5"
EGMIN="0.010"

LABEL="normtest"
NEVENTS_PER_JOB="20000"
NJOBS="50"
WALLTIME="24hr"
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
        -disk 1gb \
        -ram 256mb \
        -time "$WALLTIME" \
        -os el9 \
        -input "$INPUT" "$INPUT" \
        -- "$SCRIPT" "$INPUT" "$DIR"
done

swif2 run "$WORKFLOW"

echo "Submitted ${WORKFLOW}"
echo
echo "Check with:"
echo "  swif2 status ${WORKFLOW} -summary -problems"
echo "  swif2 status ${WORKFLOW} -jobs | head"
