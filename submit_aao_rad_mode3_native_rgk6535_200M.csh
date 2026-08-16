#!/bin/tcsh -f

# Submit globally unweighted RGK 6.535-GeV radiative e p pi0 events.
#
# Mode 3 samples (1/Q_l^2, x_l, -t_h, phi_h), performs a short envelope
# prescan in every job, and then uses AAO's stochastic multiplicity rule.
# The direct and final-LUND target boxes are the nominal RGK analysis ranges
# padded by 3.5% of the full Q2, xB, and -t spans on each side.

set script_dir = `dirname "$0"`
set script_dir = `cd "$script_dir" && pwd`
set runner = "$script_dir/run_swif_job.sh"

if ( ! -x "$runner" ) then
    echo "ERROR: missing executable SWIF runner: $runner"
    exit 1
endif

# Require an explicit campaign scale.  --full launches 40,000 SWIF jobs.
if ( $#argv != 1 || \
     ( "$argv[1]" != "--smoke" && \
       "$argv[1]" != "--pilot" && \
       "$argv[1]" != "--full" ) ) then
    echo "Usage: $0 --smoke|--pilot|--full"
    echo "  --smoke submits 1 x 100 events"
    echo "  --pilot submits 10 x 5,000 events"
    echo "  --full  submits 40,000 x 5,000 events = 200M events"
    exit 2
endif

# ------------------------- user-adjustable settings -------------------------
if ( "$argv[1]" == "--smoke" ) then
    set workflow = "aao_rad_rgk6535_mode3_native_smoke"
    set events_per_job = 100
    set jobs = 1
else if ( "$argv[1]" == "--pilot" ) then
    set workflow = "aao_rad_rgk6535_mode3_native_pilot10"
    set events_per_job = 5000
    set jobs = 10
else
    set workflow = "aao_rad_rgk6535_mode3_native_200M"
    set events_per_job = 5000
    set jobs = 40000
endif
set outbase = "/volatile/clas12/$USER/rad_mode3/$workflow"

set beam = "6.535"

# Nominal Q2 1.0:6.535, xB 0.05:0.70, -t 0.09:2.0;
# each nonperiodic range is padded by 3.5% of its full nominal span.
set q2_low = "0.806275"
set q2_high = "6.728725"
set xb_low = "0.02725"
set xb_high = "0.72275"
set t_low = "0.02315"
set t_high = "2.06685"
set phi_low = "0.0"
set phi_high = "360.0"

set w_min = "2.0"
set electron_p_min = "1.0"
set photon_e_min = "0.010"

set scan_trials = 10000
set scan_factor = "1.5"
set seed_base = 1507000001
set heartbeat = 100000

set walltime = "24hr"
set ram = "2gb"
set disk = "2gb"
set files_per_directory = 5000
# ---------------------------------------------------------------------------

@ total_events = $events_per_job * $jobs
if ( "$argv[1]" == "--full" && $total_events != 200000000 ) then
    echo "ERROR: events_per_job * jobs must equal 200000000"
    exit 2
endif

mkdir -p "$outbase/inputs"

echo "Creating $workflow"
echo "  total events: $total_events ($jobs jobs x $events_per_job)"
echo "  direct proposal: 1/Q_l^2, x_l, -t_h, phi_h"
echo "  padded target: Q2=${q2_low}:${q2_high}, xB=${xb_low}:${xb_high}"
echo "                 -t=${t_low}:${t_high}, phi=${phi_low}:${phi_high}"
echo "  W > $w_min GeV, final electron p > $electron_p_min GeV"
echo "  no separate y_max cut"
echo "  automatic scan: $scan_trials trials/job, factor=$scan_factor"
echo "  output: $outbase"

swif2 create -workflow "$workflow"

@ task = 0
while ( $task < $jobs )
    set task_tag = `printf "%08d" $task`
    @ chunk_index = $task / $files_per_directory
    set chunk_tag = `printf "%04d" $chunk_index`
    @ seed = $seed_base + $task
    set input_dir = "$outbase/inputs/chunk_$chunk_tag"
    set output_dir = "$outbase/lund/chunk_$chunk_tag"
    set input = "$input_dir/mode3_rgk6535__g${task_tag}.inp"

    mkdir -p "$input_dir" "$output_dir"

    cat >! "$input" << EOF
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
$beam
$q2_low $q2_high
$electron_p_min $beam
$photon_e_min
$events_per_job
$scan_factor
3
$seed
$task
1.0
$xb_low $xb_high
$t_low $t_high
$phi_low $phi_high
$q2_low $q2_high
$xb_low $xb_high
$t_low $t_high
$phi_low $phi_high
$w_min
0 1.0
2
$scan_trials
$heartbeat
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

    @ task++
end


swif2 run "$workflow"

echo "Submitted $workflow"
echo "Monitor with:"
echo "  swif2 diagnose $workflow"
echo "  swif2 status $workflow"
echo "  swif2 status $workflow --problems"

