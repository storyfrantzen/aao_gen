#!/bin/bash

# ---- strict mode ----
set -euo pipefail

echo "==== SWIF JOB START ===="

# ---- job identity ----
UNIQ_ID="${SWIF_JOB_ID:-manual_$$}"
echo "Job ID: ${UNIQ_ID}"
echo "User: ${USER}"
echo "Host: $(hostname)"
echo "Start time: $(date)"

# ---- modules ----
if [ -f /etc/profile.d/modules.sh ]; then
    source /etc/profile.d/modules.sh
else
    echo "ERROR: modules system not found"
    exit 127
fi

module use /cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/modulefiles
module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
module load clas12/5.4

echo "Modules loaded"
module list

# ---- args ----
if [ $# -lt 2 ]; then
    echo "ERROR: Usage: run_swif_job.sh <input.inp> <output_dir>"
    exit 1
fi

INPUT_ABS="$1"
OUTPUT_BASE="$2"

AAO_EXE="/work/clas12/storyf/aao_gen/aao_rad/build/aao_rad.exe"

echo "Executable: $AAO_EXE"
echo "Input file: $INPUT_ABS"
echo "Output base: $OUTPUT_BASE"

# ---- checks ----
if [ ! -x "$AAO_EXE" ]; then
    echo "ERROR: executable not found: $AAO_EXE"
    exit 127
fi

if [ ! -f "$INPUT_ABS" ]; then
    echo "ERROR: input not found: $INPUT_ABS"
    exit 1
fi

# ---- extract kinematics (robust positional parsing) ----
# (based on your current input layout)

EBEAM=$(sed -n '12p' "$INPUT_ABS" | awk '{print $1}')
Q2MIN=$(sed -n '13p' "$INPUT_ABS" | awk '{print $1}')
Q2MAX=$(sed -n '13p' "$INPUT_ABS" | awk '{print $2}')
EPMIN=$(sed -n '14p' "$INPUT_ABS" | awk '{print $1}')
EPMAX=$(sed -n '14p' "$INPUT_ABS" | awk '{print $2}')
EGMIN=$(sed -n '15p' "$INPUT_ABS" | awk '{print $1}')

# ---- physics tag ----
TAG="E${EBEAM}_Q2${Q2MIN}-${Q2MAX}_EP${EPMIN}-${EPMAX}_EG${EGMIN}"

echo "Kinematic tag: $TAG"

# ---- workspace ----
WORKROOT="${SWIF_JOB_WORK_DIR:-/tmp}"
WORKDIR="${WORKROOT}/aao_${UNIQ_ID}"
mkdir -p "$WORKDIR"
cleanup() {
    cd / || true
    rm -rf "$WORKDIR"
}
trap cleanup EXIT
cd "$WORKDIR"

echo "Working directory: $WORKDIR"

# ---- stage inputs ----
cp "$AAO_EXE" ./aao_rad.exe
cp "$INPUT_ABS" ./input.inp

echo "Local files:"
ls -lh

# ---- run generator ----
echo "Starting aao_rad.exe"
set +e
./aao_rad.exe < input.inp
STATUS=$?
set -e

if [ $STATUS -ne 0 ]; then
    echo "ERROR: generator failed with code $STATUS"
    exit 2
fi

# ---- output check ----
if [ ! -f aao_rad.lund ]; then
    echo "ERROR: missing output file aao_rad.lund"
    exit 3
fi

# ---- output destination ----
mkdir -p "$OUTPUT_BASE"
OUTSTEM="${OUTPUT_BASE}/aao_rad_${TAG}_${UNIQ_ID}"

mv aao_rad.lund "${OUTSTEM}.lund"

if [ -f aao_rad.norm ]; then
    mv aao_rad.norm "${OUTSTEM}.norm"
else
    echo "WARNING: missing output file aao_rad.norm"
fi

if [ -f aao_rad.sum ]; then
    mv aao_rad.sum "${OUTSTEM}.sum"
else
    echo "WARNING: missing output file aao_rad.sum"
fi

if [ -f aao_rad.out ]; then
    mv aao_rad.out "${OUTSTEM}.out"
else
    echo "WARNING: missing output file aao_rad.out"
fi

cp input.inp "${OUTSTEM}.inp"

echo "Saved output stem: ${OUTSTEM}"
ls -lh "${OUTSTEM}".*

# ---- cleanup ----
cd /
rm -rf "$WORKDIR"

echo "==== SUCCESS ${UNIQ_ID} ===="
echo "End time: $(date)"
exit 0