#!/bin/tcsh -f

# Submit globally unweighted Born e p pi0 events for the padded RGK region.
#
# Production is deliberately split into eight independent 25M-event batches:
#   5,000 jobs/batch x 5,000 events/job x 8 batches = 200M events.
#
# Run one batch at a time with --batch 0 through --batch 7.  Each job requests
# 5,000 events to align with the OSG type-2 limit. AAO stochastic multiplicity
# can rarely overshoot nmax, so exact per-file counts must still be audited
# before OSG submission.

set script_dir = `dirname "$0"`
set script_dir = `cd "$script_dir" && pwd`
set runner = "$script_dir/run_swif_norad_job.sh"

if ( ! -x "$runner" ) then
    echo "ERROR: missing executable SWIF runner: $runner"
    exit 1
endif

if ( $#argv < 1 || $#argv > 2 ) then
    echo "Usage: $0 --smoke|--pilot|--batch BATCH_INDEX"
    echo "       $0 --nominal-smoke|--nominal-pilot|--nominal-batch BATCH_INDEX"
    echo "  --smoke   submits 1 x 100 events"
    echo "  --pilot   submits 10 x 5,000 events"
    echo "  --batch N submits production batch N, where N is 0 through 7"
    exit 2
endif

set mode = "$argv[1]"
set events_per_job = 5000
set first_task = 0
set jobs = 0
set production = 0
set nominal_conditioned = 0

if ( "$mode" == "--smoke" && $#argv == 1 ) then
    set workflow = "aao_norad_rgk6535_mode3_smoke"
    set events_per_job = 100
    set jobs = 1
    set campaign_root = "/volatile/clas12/$USER/norad_mode3/$workflow"
    set outbase = "$campaign_root"
else if ( "$mode" == "--pilot" && $#argv == 1 ) then
    set workflow = "aao_norad_rgk6535_mode3_pilot10"
    set jobs = 10
    set campaign_root = "/volatile/clas12/$USER/norad_mode3/$workflow"
    set outbase = "$campaign_root"
else if ( "$mode" == "--batch" && $#argv == 2 ) then
    set batch = "$argv[2]"
    if ( "$batch" !~ [0-7] ) then
        echo "ERROR: BATCH_INDEX must be one integer from 0 through 7"
        exit 2
    endif
    set production = 1
    set jobs = 5000
    @ first_task = $batch * $jobs
    set batch_tag = `printf "%02d" $batch`
    set workflow = "aao_norad_rgk6535_mode3_200M_b${batch_tag}"
    set campaign_root = "/volatile/clas12/$USER/norad_mode3/aao_norad_rgk6535_mode3_200M"
    set outbase = "$campaign_root/batch_$batch_tag"
else if ( "$mode" == "--nominal-smoke" && $#argv == 1 ) then
    set nominal_conditioned = 1
    set workflow = "aao_norad_rgk6535_mode3_nominal_w2_smoke"
    set events_per_job = 100
    set jobs = 1
    set campaign_root = "/volatile/clas12/$USER/norad_mode3/$workflow"
    set outbase = "$campaign_root"
else if ( "$mode" == "--nominal-pilot" && $#argv == 1 ) then
    set nominal_conditioned = 1
    set workflow = "aao_norad_rgk6535_mode3_nominal_w2_pilot10"
    set jobs = 10
    set campaign_root = "/volatile/clas12/$USER/norad_mode3/$workflow"
    set outbase = "$campaign_root"
else if ( "$mode" == "--nominal-batch" && $#argv == 2 ) then
    set batch = "$argv[2]"
    if ( "$batch" !~ [0-7] ) then
        echo "ERROR: BATCH_INDEX must be one integer from 0 through 7"
        exit 2
    endif
    set nominal_conditioned = 1
    set production = 1
    set jobs = 5000
    @ first_task = $batch * $jobs
    set batch_tag = `printf "%02d" $batch`
    set workflow = "aao_norad_rgk6535_mode3_nominal_w2_200M_b${batch_tag}"
    set campaign_root = "/volatile/clas12/$USER/norad_mode3/aao_norad_rgk6535_mode3_nominal_w2_200M"
    set outbase = "$campaign_root/batch_$batch_tag"
else
    echo "Usage: $0 --smoke|--pilot|--batch BATCH_INDEX"
    echo "       $0 --nominal-smoke|--nominal-pilot|--nominal-batch BATCH_INDEX"
    exit 2
endif

# RGK analysis binning padded by 3.5% of each complete nonperiodic span.
set beam = "6.535"
set q2_low = "0.806275"
set q2_high = "6.728725"
set xb_low = "0.02725"
set xb_high = "0.72275"
set t_low = "0.02315"
set t_high = "2.06685"
set phi_low = "0.0"
set phi_high = "360.0"
set electron_p_min = "1.0"
set w_min = "2.0"
set condition_phase_space = 0
if ( $nominal_conditioned ) then
    set q2_low = "1.0"
    set q2_high = "6.535"
    set xb_low = "0.05"
    set xb_high = "0.70"
    set t_low = "0.09"
    set t_high = "2.0"
    set condition_phase_space = 1
endif

set physics_model = 5
set fmcall = "2.0"
set seed_base = 1807000001

# Conservative requests based on the generator's small observed memory use.
set walltime = "8hr"
set ram = "256mb"
set disk = "1gb"
set files_per_directory = 2500

@ total_events = $events_per_job * $jobs
@ last_task = $first_task + $jobs - 1

mkdir -p "$outbase/inputs" "$outbase/lund"

echo "Creating $workflow"
echo "  task range: $first_task through $last_task"
echo "  events: $total_events ($jobs jobs x $events_per_job)"
echo "  proposal: 1/Q2, xB, -t, phi (Born mode 3)"
if ( $nominal_conditioned ) then
    echo "  nominal conditioned box: Q2=${q2_low}:${q2_high}, xB=${xb_low}:${xb_high}"
else
    echo "  padded box: Q2=${q2_low}:${q2_high}, xB=${xb_low}:${xb_high}"
endif
echo "              -t=${t_low}:${t_high}, phi=${phi_low}:${phi_high}"
echo "  electron momentum: ${electron_p_min}:${beam} GeV"
echo "  phase-space conditioning flag: $condition_phase_space"
if ( $condition_phase_space ) echo "  generated events require W >= $w_min GeV; no y_max below 1"
echo "  resources/job: 1 core, $ram RAM, $disk disk, $walltime"
echo "  output: $outbase"
if ( ! $condition_phase_space ) then
    echo "  NOTE: apply the nominal W >= 2 GeV analysis cut downstream."
endif

swif2 create -workflow "$workflow"

@ local_task = 0
while ( $local_task < $jobs )
    @ global_task = $first_task + $local_task
    @ chunk_index = $local_task / $files_per_directory
    @ seed_positive = $seed_base + $global_task
    @ seed = -1 * $seed_positive

    set task_tag = `printf "%08d" $global_task`
    set chunk_tag = `printf "%04d" $chunk_index`
    set input_dir = "$outbase/inputs/chunk_$chunk_tag"
    set output_dir = "$outbase/lund/chunk_$chunk_tag"
    set input = "$input_dir/mode3_born_rgk6535__g${task_tag}.inp"

    mkdir -p "$input_dir" "$output_dir"

    cat >! "$input" << EOF
$physics_model
1
3
1
$beam
$q2_low $q2_high
$electron_p_min $beam
$events_per_job
$fmcall
0
$seed
3
$xb_low $xb_high
$t_low $t_high
$phi_low $phi_high
$condition_phase_space $w_min 1.0
EOF

    swif2 add-job \
        -workflow "$workflow" \
        -name "${workflow}_g${task_tag}" \
        -cores 1 \
        -ram "$ram" \
        -disk "$disk" \
        -time "$walltime" \
        -os el9 \
        -input "$input" "$input" \
        -- "$runner" "$input" "$output_dir"

    @ local_task++
end

swif2 run "$workflow"

echo "Submitted $workflow"
if ( $production ) then
    echo "This is production batch $batch_tag of 00 through 07."
endif
echo "Monitor with:"
echo "  swif2 diagnose $workflow"
echo "  swif2 status $workflow"
echo "  swif2 status $workflow --problems"
