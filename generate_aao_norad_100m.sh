#!/bin/bash

# Generate the RGK 6.535 GeV AAO Born/non-radiative sample in local chunks.
#
# Default target:
#   100M generated events = 20000 chunks * 5000 events/chunk
#
# Run inside tmux on ifarm:
#   cd /work/clas12/storyf/SF_analysis_software_v2.0/external/aao_print/aao_gen
#   tmux new -s aao_norad_100M
#   ./generate_aao_norad_100m.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AAO_GEN_DIR="${AAO_GEN_DIR:-${SCRIPT_DIR}}"
AAO_NORAD_DIR="${AAO_NORAD_DIR:-${AAO_GEN_DIR}/aao_norad}"
AAO_EXE="${AAO_EXE:-${AAO_NORAD_DIR}/build/aao_norad}"

OUTDIR="${OUTDIR:-/volatile/clas12/${USER}/LUND/aao_norad_100M/aao_norad_E6.535_Q2_0.9_5.0_EP_1.1_5.0}"
WORKROOT="${SWIF_JOB_WORK_DIR:-${TMPDIR:-/scratch/${USER}/tmpfs}}"

BEAM="${BEAM:-6.535}"
Q2MIN="${Q2MIN:-0.9}"
Q2MAX="${Q2MAX:-5.0}"
EPMIN="${EPMIN:-1.1}"
EPMAX="${EPMAX:-5.0}"
EVENTS_PER_CHUNK="${EVENTS_PER_CHUNK:-5000}"
NCHUNKS="${NCHUNKS:-20000}"
FIRST_CHUNK="${FIRST_CHUNK:-1}"

PHYS="${PHYS:-5}"
FLAG_EHEL="${FLAG_EHEL:-1}"
NPART="${NPART:-3}"
EPIREA="${EPIREA:-1}"
FMCALL="${FMCALL:-1.0}"
BOSO="${BOSO:-1}"
SEED_SOURCE="${SEED_SOURCE:-0}"

echo "==== AAO NORAD LOCAL GENERATION START ===="
echo "AAO generator dir: ${AAO_GEN_DIR}"
echo "Executable: ${AAO_EXE}"
echo "Output directory: ${OUTDIR}"
echo "Work root: ${WORKROOT}"
echo "Chunks: ${FIRST_CHUNK}-${NCHUNKS}"
echo "Events/chunk: ${EVENTS_PER_CHUNK}"
echo "Total requested events from chunk range: $(( (NCHUNKS - FIRST_CHUNK + 1) * EVENTS_PER_CHUNK ))"
echo "Kinematics: E=${BEAM}, Q2=${Q2MIN}-${Q2MAX}, Ep=${EPMIN}-${EPMAX}"
echo "Start time: $(date)"

mkdir -p "${OUTDIR}" "${WORKROOT}"

if [ ! -x "${AAO_EXE}" ]; then
    echo "Executable not found, attempting build in ${AAO_NORAD_DIR}"
    make -C "${AAO_NORAD_DIR}"
fi

if [ ! -x "${AAO_EXE}" ]; then
    echo "ERROR: executable not found or not executable: ${AAO_EXE}" >&2
    exit 127
fi

cleanup_workdir() {
    local workdir="$1"
    if [ -n "${workdir}" ] && [ -d "${workdir}" ]; then
        rm -rf "${workdir}"
    fi
}

for job in $(seq "${FIRST_CHUNK}" "${NCHUNKS}"); do
    JOB_ID="$(printf "%05d" "${job}")"
    WORKDIR="${WORKROOT}/aao_norad_${JOB_ID}_$$"
    OUTSTEM="${OUTDIR}/aao_norad_E${BEAM}_Q2${Q2MIN}-${Q2MAX}_EP${EPMIN}-${EPMAX}_${JOB_ID}"

    if [ -f "${OUTSTEM}.norm" ] && [ -f "${OUTSTEM}.lund" ]; then
        echo "Skipping ${JOB_ID}; ${OUTSTEM}.norm and .lund already exist"
        continue
    fi

    cleanup_workdir "${WORKDIR}"
    mkdir -p "${WORKDIR}"

    cat > "${WORKDIR}/input.inp" <<EOF
${PHYS}
${FLAG_EHEL}
${NPART}
${EPIREA}
${BEAM}
${Q2MIN} ${Q2MAX}
${EPMIN} ${EPMAX}
${EVENTS_PER_CHUNK}
${FMCALL}
${BOSO}
${SEED_SOURCE}
EOF

    cp "${AAO_EXE}" "${WORKDIR}/aao_norad.exe"

    (
        cd "${WORKDIR}"
        ./aao_norad.exe < input.inp > aao_norad.stdout 2>&1

        test -f aao_norad.lund
        test -f aao_norad.norm

        mv aao_norad.lund "${OUTSTEM}.lund"
        mv aao_norad.norm "${OUTSTEM}.norm"
        [ -f aao_norad.kin ] && mv aao_norad.kin "${OUTSTEM}.kin"
        [ -f aao_norad.sum ] && mv aao_norad.sum "${OUTSTEM}.sum"
        [ -f aao_norad.out ] && mv aao_norad.out "${OUTSTEM}.out"
        cp input.inp "${OUTSTEM}.inp"
        mv aao_norad.stdout "${OUTSTEM}.stdout"
    )

    cleanup_workdir "${WORKDIR}"

    if (( job % 100 == 0 )); then
        echo "$(date): completed chunk ${job}/${NCHUNKS}"
    fi
done

echo "==== AAO NORAD LOCAL GENERATION DONE ===="
echo "End time: $(date)"
echo "Output directory: ${OUTDIR}"
echo "LUND files:"
find "${OUTDIR}" -maxdepth 1 -type f -name '*.lund' | wc -l
echo "Norm files:"
find "${OUTDIR}" -maxdepth 1 -type f -name '*.norm' | wc -l
