#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="${SCRIPT_DIR}/run_swif_norad_job.sh"

OUTBASE="/volatile/clas12/storyf/LUND/aao_norad_normtest"

BEAM="6.535"
Q2MIN="0.9"
Q2MAX="5.0"
EPMIN="1.1"
EPMAX="5.0"

LABEL="normtest"
NEVENTS_PER_JOB="5000"
NJOBS="200"
WALLTIME="24hr"

WORKFLOW="aao_norad_${LABEL}_E${BEAM}_Q2_${Q2MIN}_${Q2MAX}_EP_${EPMIN}_${EPMAX}"
DIR="${OUTBASE}/${WORKFLOW}"
INPUT_DIR="${DIR}/inputs"

mkdir -p "$INPUT_DIR"

swif2 create -workflow "$WORKFLOW"

for job in $(seq 1 "$NJOBS"); do
    JOB_ID=$(printf "%05d" "$job")
    INPUT="${INPUT_DIR}/input_${JOB_ID}.inp"

    cat > "$INPUT" << EOF_INP
5
1
3
1
${BEAM}
${Q2MIN} ${Q2MAX}
${EPMIN} ${EPMAX}
${NEVENTS_PER_JOB}
1.0
1
0
EOF_INP

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
echo "Check with:"
echo "  swif2 status ${WORKFLOW} -summary -problems"
