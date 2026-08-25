#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: $0 BATCH_INDEX [PARALLEL_PROCESSES]" >&2
    echo "  BATCH_INDEX: 0 through 7 (25M requested events per batch)" >&2
    echo "  PARALLEL_PROCESSES: defaults to 4" >&2
}

if [[ $# -lt 1 || $# -gt 2 || ! $1 =~ ^[0-7]$ ]]; then
    usage
    exit 2
fi

BATCH=$1
PARALLEL=${2:-4}
if [[ ! $PARALLEL =~ ^[1-9][0-9]*$ || $PARALLEL -gt 16 ]]; then
    echo "ERROR: PARALLEL_PROCESSES must be an integer from 1 through 16" >&2
    exit 2
fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
RUNNER="$SCRIPT_DIR/run_swif_norad_job.sh"
AAO_EXE="$SCRIPT_DIR/aao_norad/build/aao_norad"
CAMPAIGN=${AAO_NORAD_MANUAL_CAMPAIGN:-/volatile/clas12/$USER/norad_mode3/aao_norad_rgk6535_mode3_nominal_w2_200M}
SCRATCH=${AAO_NORAD_MANUAL_SCRATCH:-/scratch/$USER/tmpfs}

if [[ ! -x $RUNNER || ! -x $AAO_EXE ]]; then
    echo "ERROR: build aao_norad and verify the runner before starting:" >&2
    echo "  make -C $SCRIPT_DIR/aao_norad -j2" >&2
    exit 3
fi

printf -v BATCH_TAG '%02d' "$BATCH"
BATCH_DIR="$CAMPAIGN/batch_$BATCH_TAG"
TASK_FILE="$BATCH_DIR/task_indices.txt"
mkdir -p "$BATCH_DIR/inputs" "$BATCH_DIR/lund" "$BATCH_DIR/logs" "$SCRATCH"

exec 9>"$BATCH_DIR/.manual_batch.lock"
if ! flock -n 9; then
    echo "ERROR: another manual runner already holds the batch-$BATCH_TAG lock" >&2
    exit 4
fi

exec > >(tee -a "$BATCH_DIR/manual_batch.log") 2>&1

EVENTS_PER_JOB=5000
JOBS=5000
FILES_PER_DIRECTORY=2500
SEED_BASE=1907000001
FIRST_TASK=$((BATCH * JOBS))

echo "[$(date -Is)] Preparing RGK Born nominal-W2 batch $BATCH_TAG"
echo "  parallel processes: $PARALLEL"
echo "  task range: $FIRST_TASK through $((FIRST_TASK + JOBS - 1))"
echo "  output: $BATCH_DIR"

: > "$TASK_FILE"
for ((local_task = 0; local_task < JOBS; local_task++)); do
    global_task=$((FIRST_TASK + local_task))
    chunk_index=$((local_task / FILES_PER_DIRECTORY))
    seed=$((-(SEED_BASE + global_task)))
    printf -v task_tag '%08d' "$global_task"
    printf -v chunk_tag '%04d' "$chunk_index"
    input_dir="$BATCH_DIR/inputs/chunk_$chunk_tag"
    output_dir="$BATCH_DIR/lund/chunk_$chunk_tag"
    input="$input_dir/mode3_born_rgk6535__g${task_tag}.inp"
    mkdir -p "$input_dir" "$output_dir"
    cat > "$input" <<EOF
5
1
3
1
6.535
1.0 6.535
1.0 6.535
$EVENTS_PER_JOB
2.0
0
$seed
3
0.05 0.70
0.09 2.0
0.0 360.0
1 2.0 1.0
EOF
    echo "$global_task" >> "$TASK_FILE"
done

export SCRIPT_DIR CAMPAIGN SCRATCH BATCH_DIR
run_one() {
    local global_task=$1
    local batch=$((global_task / 5000))
    local local_task=$((global_task % 5000))
    local chunk=$((local_task / 2500))
    local batch_tag task_tag chunk_tag input output_dir log stem
    printf -v batch_tag '%02d' "$batch"
    printf -v task_tag '%08d' "$global_task"
    printf -v chunk_tag '%04d' "$chunk"
    input="$CAMPAIGN/batch_$batch_tag/inputs/chunk_$chunk_tag/mode3_born_rgk6535__g${task_tag}.inp"
    output_dir="$CAMPAIGN/batch_$batch_tag/lund/chunk_$chunk_tag"
    log="$CAMPAIGN/batch_$batch_tag/logs/task_${task_tag}.log"
    stem="$output_dir/aao_norad_E6.535_Q21.0-6.535_EP1.0-6.535_g${task_tag}"

    if [[ -s ${stem}.lund && -s ${stem}.norm ]]; then
        echo "SKIP completed task $task_tag"
        return 0
    fi

    echo "START task $task_tag"
    if TMPDIR="$SCRATCH" SWIF_JOB_ID="g${task_tag}" \
        "$SCRIPT_DIR/run_swif_norad_job.sh" "$input" "$output_dir" \
        > "$log" 2>&1; then
        echo "DONE task $task_tag"
    else
        local status=$?
        echo "FAILED task $task_tag (status $status); inspect $log"
        return "$status"
    fi
}
export -f run_one

set +e
xargs -P "$PARALLEL" -n 1 bash -c 'run_one "$1"' _ < "$TASK_FILE"
status=$?
set -e

norm_count=$(find "$BATCH_DIR/lund" -type f -name '*.norm' | wc -l)
lund_count=$(find "$BATCH_DIR/lund" -type f -name '*.lund' | wc -l)
echo "[$(date -Is)] Batch $BATCH_TAG stopped with status $status"
echo "  completed norm files: $norm_count / $JOBS"
echo "  completed LUND files: $lund_count / $JOBS"
exit "$status"
