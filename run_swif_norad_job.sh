#!/bin/bash
set -euo pipefail

echo "==== NORAD SWIF JOB START ===="
UNIQ_ID="${SWIF_JOB_ID:-manual_$$}"

if [ -f /etc/profile.d/modules.sh ]; then
    source /etc/profile.d/modules.sh
else
    echo "ERROR: modules system not found"
    exit 127
fi

module use /cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/modulefiles
module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
module load clas12/5.4

if [ $# -lt 2 ]; then
    echo "Usage: run_swif_norad_job.sh <input.inp> <output_dir>"
    exit 1
fi

INPUT_ABS="$1"
OUTPUT_BASE="$2"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AAO_EXE="${SCRIPT_DIR}/aao_norad/build/aao_norad"

if [ ! -x "$AAO_EXE" ]; then
    echo "ERROR: executable not found: $AAO_EXE"
    exit 127
fi

EBEAM=$(sed -n '5p' "$INPUT_ABS" | awk '{print $1}')
Q2MIN=$(sed -n '6p' "$INPUT_ABS" | awk '{print $1}')
Q2MAX=$(sed -n '6p' "$INPUT_ABS" | awk '{print $2}')
EPMIN=$(sed -n '7p' "$INPUT_ABS" | awk '{print $1}')
EPMAX=$(sed -n '7p' "$INPUT_ABS" | awk '{print $2}')

TAG="E${EBEAM}_Q2${Q2MIN}-${Q2MAX}_EP${EPMIN}-${EPMAX}"

WORKROOT="${SWIF_JOB_WORK_DIR:-${TMPDIR:-/tmp}}"
WORKDIR="${WORKROOT}/aao_norad_${UNIQ_ID}"
mkdir -p "$WORKDIR"

cleanup() {
    cd / || true
    rm -rf "$WORKDIR"
}
trap cleanup EXIT

cd "$WORKDIR"
cp "$AAO_EXE" ./aao_norad.exe
cp "$INPUT_ABS" ./input.inp

./aao_norad.exe < input.inp

if [ ! -f aao_norad.lund ]; then
    echo "ERROR: missing aao_norad.lund"
    exit 3
fi

mkdir -p "$OUTPUT_BASE"
OUTSTEM="${OUTPUT_BASE}/aao_norad_${TAG}_${UNIQ_ID}"

mv aao_norad.lund "${OUTSTEM}.lund"

[ -f aao_norad.norm ] && mv aao_norad.norm "${OUTSTEM}.norm"
[ -f aao_norad.sum ]  && mv aao_norad.sum  "${OUTSTEM}.sum"
[ -f aao_norad.out ]  && mv aao_norad.out  "${OUTSTEM}.out"

cp input.inp "${OUTSTEM}.inp"

echo "Saved output stem: ${OUTSTEM}"
ls -lh "${OUTSTEM}".*
echo "==== NORAD SUCCESS ${UNIQ_ID} ===="
