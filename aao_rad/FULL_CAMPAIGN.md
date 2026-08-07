# Full radiative mode-4 campaign

`radiative_full_campaign.py` promotes the tested per-stratum mode-4 workflow
to an audited, scheduler-neutral campaign.  It does not replace any physics
or normalization calculation in `radiative_mode4.py`.  Instead, it prepares
the existing immutable manifests, assigns one run per farm task, verifies
completed artifacts, pools calibration campaigns, and sizes follow-ups from
the measured target rates.

The default RGK campaign operates on the selected cumulative queue.  For the
current `y_max = 0.95` census this means all 2,301 data-occupied strata plus
the 1,336 zero-data model strata needed to leave only a 0.1% global model
residual: 3,637 strata total.  The formal 12,960-bin Cartesian catalog can be
requested with `--selection all-catalog`, but that is not the production
default because it includes bins with no demonstrated support.

## Safety and reproducibility

- Planning snapshots the analysis configuration, continuous-guard recipes,
  legacy generator input, optional refinements, and flat-index selection.
- Every task freezes its manifest SHA-256, flat index, replica, seed, output
  location, and requested exposure.
- `run-task` is resumable.  A valid completed task is reported as
  `already_complete`; a failed task writes a traceback-bearing receipt and
  can be retried by SWIF.
- `status` counts only run JSON files whose schema, stratum, replica, and
  source-manifest hash agree with the task table.
- `emit-swif` writes a submission script but never invokes `swif2` itself.
- Follow-up calibration manifests remain compatible with the initial
  manifest and are pooled component by component by the existing mode-4
  finalizer.
- Production planning fails closed unless every selected queue stratum has
  a usable envelope.  `--allow-incomplete` is available only as an explicit
  development override.

## 1. Plan the complete RGK calibration

From `aao_rad` on the farm:

```tcsh
set repo = /w/hallb-scshelf2102/clas12/storyf/SF_analysis_software_v2.0
set aao_rad = "$repo/external/aao_gen/aao_rad"
set fresh = "$aao_rad/rgk_mode4_v3_fresh_20260730"
set census = "$fresh/census/rgk_iteration000"

set census_config = "$census/inputs/analysis_ymax0p95.json"
set queue = "$census/outputs/cumulative_queue_iteration000"
set continuous = "$census/outputs/continuous_guards_iteration000"
set legacy_input = "$fresh/batches/representative_batch_001/inputs/legacy_input.inp"
set campaign = "$census/campaigns/rgk_full_calibration_iteration000"

cd "$aao_rad"
set revision = `git rev-parse HEAD`
```

If the frozen configuration has a different filename, point
`census_config` at the file whose SHA-256 is
`83c3d905fc4780320ad67dae2dbde06e75de0eaa24386ff51cfdc8d2367504a2`.

```tcsh
python3 radiative_full_campaign.py plan-calibration \
  --queue "$queue/cumulative_stratum_queue.json" \
  --config "$census_config" \
  --recipes "$continuous/continuous_guard_recipes.json" \
  --input "$legacy_input" \
  --output "$campaign" \
  --selection selected \
  --candidate padding_0p035 \
  --core-fraction 0.90 \
  --inside-guard-trial-fraction 0.50 \
  --trials 7000000 \
  --replicas 1 \
  --heartbeat-interval 100000 \
  --apply-y-max \
  --seed-base 907001 \
  --generator-revision "$revision"
```

Verify the plan:

```tcsh
python3 -c 'import json,sys; c=json.load(open(sys.argv[1])); print("kind =",c["kind"]); print("tasks =",c["task_count"]); print("selection =",json.dumps(c["selection"],indent=2))' \
  "$campaign/campaign.json"

wc -l "$campaign/task_ids.txt"
```

The current selected queue should produce 3,637 tasks.

## 2. Emit and submit the SWIF2 workflow

```tcsh
python3 radiative_full_campaign.py emit-swif \
  --campaign "$campaign/campaign.json" \
  --workflow rgk_mode4_y095_calibration_i000 \
  --executable "$aao_rad/build/aao_rad" \
  --walltime 8hr \
  --ram 1gb \
  --disk 2gb
```

Inspect before submission:

```tcsh
less "$campaign/submit_swif.sh"
```

Submit only after inspection:

```tcsh
"$campaign/submit_swif.sh"
```

SWIF monitoring remains external to the driver:

```tcsh
swif2 status rgk_mode4_y095_calibration_i000 -summary -problems
swif2 status rgk_mode4_y095_calibration_i000 -jobs | head
```

The emitted workflow is safely resumable.  Resubmitting a task that already
has a valid run artifact returns `already_complete` without rerunning AAO.

## 3. Collect status and finalize

```tcsh
python3 radiative_full_campaign.py status \
  --campaign "$campaign/campaign.json"

python3 -c 'import json,sys; s=json.load(open(sys.argv[1])); print(json.dumps(s["status_counts"],indent=2)); print("complete =",s["complete"])' \
  "$campaign/status.json"
```

When all tasks are complete:

```tcsh
python3 radiative_full_campaign.py finalize \
  --campaign "$campaign/campaign.json" \
  --envelope-safety-factor 1.20 \
  --maximum-duplicate-fraction 0.05 \
  --minimum-component-targets 20 \
  --minimum-provisional-inside-targets 1000 \
  --allow-zero-complement \
  --zero-complement-confidence 0.95 \
  --maximum-zero-complement-target-rate 1e-6
```

This writes `envelope_calibration.json` and
`envelope_calibration.tsv` under the campaign root.

## 4. Automatically prepare follow-ups

Count readiness:

```tcsh
python3 -c 'import csv,collections,sys; r=list(csv.DictReader(open(sys.argv[1]),delimiter="\t")); print("strata =",len(r)); print("readiness =",dict(collections.Counter(x["pilot_readiness"] for x in r))); print("status =",dict(collections.Counter(x["recommendation_status"] for x in r)))' \
  "$campaign/envelope_calibration.tsv"
```

If any strata are not ready:

```tcsh
set followup = "$census/campaigns/rgk_full_calibration_followup_iteration001"

python3 radiative_full_campaign.py plan-followup \
  --campaign "$campaign/campaign.json" \
  --calibration "$campaign/envelope_calibration.json" \
  --output "$followup" \
  --minimum-component-targets 20 \
  --minimum-provisional-inside-targets 1000 \
  --zero-complement-confidence 0.95 \
  --maximum-zero-complement-target-rate 1e-6 \
  --followup-target-safety-factor 1.5 \
  --trial-quantum 1000000 \
  --maximum-followup-trials 100000000 \
  --discovery-trials 20000000
```

The planner groups strata by need and exposure:

- `inside_support`: 90% inside-guard trials;
- `complement_support`: 90% guard-complement trials;
- `zero_complement_exposure`: enough complement trials to reach the exact
  zero-success occurrence-rate criterion;
- `both_component_support`: 50/50 exposure;
- `discovery`: a bounded 50/50 campaign when neither component has support;
- `manual_envelope_review`: recorded but not silently scheduled when more
  calibration cannot resolve an envelope/duplicate conflict.

Inspect `followup/followup_plan.tsv`, then emit and submit its SWIF workflow:

```tcsh
python3 radiative_full_campaign.py emit-swif \
  --campaign "$followup/campaign.json" \
  --workflow rgk_mode4_y095_calibration_i001 \
  --executable "$aao_rad/build/aao_rad" \
  --walltime 24hr

less "$followup/submit_swif.sh"
"$followup/submit_swif.sh"
```

After completion, finalize the follow-up campaign.  Its campaign metadata
already includes both the initial and follow-up manifests, so the new report
is pooled automatically:

```tcsh
python3 radiative_full_campaign.py finalize \
  --campaign "$followup/campaign.json" \
  --allow-zero-complement
```

Repeat `plan-followup` only if the pooled report still contains non-ready
strata.  No manual list construction is required.

## 5. Plan generated LUND pilots or production

Once every selected stratum is ready, prepare two independent 200-event
pilots per stratum:

Use the latest follow-up campaign and its pooled report if follow-ups were
needed.  If the initial report was already complete, use `campaign` and its
report instead.

```tcsh
set final_campaign = "$followup"
set final_envelopes = "$followup/envelope_calibration.json"
set pilots = "$census/campaigns/rgk_full_pilots_2x200_iteration000"

python3 radiative_full_campaign.py plan-production \
  --campaign "$final_campaign/campaign.json" \
  --calibration "$final_envelopes" \
  --output "$pilots" \
  --events-per-stratum 200 \
  --replicas 2 \
  --seed-base 307000001
```

Emit and submit exactly as for calibration, using a distinct workflow name
and a longer walltime for sparse strata:

```tcsh
python3 radiative_full_campaign.py emit-swif \
  --campaign "$pilots/campaign.json" \
  --workflow rgk_mode4_y095_pilots_2x200_i000 \
  --executable "$aao_rad/build/aao_rad" \
  --walltime 24hr
```

After all pilot tasks complete:

```tcsh
python3 radiative_full_campaign.py validate-pilots \
  --campaign "$pilots/campaign.json" \
  --calibration "$final_envelopes" \
  --minimum-runs 2 \
  --minimum-events 400 \
  --maximum-duplicate-fraction 0.05 \
  --maximum-guard-complement-fraction 0.02 \
  --maximum-relative-cross-section-difference 0.10 \
  --maximum-cross-section-z-score 3
```

For type-2 OSG production, prepare at most 5,000 events per stratum/replica.
The resulting files retain the existing structure and naming convention:

```text
stages/production/runs/s04608/s04608__g0000.lund
```

Use additional replicas for more than 5,000 events per stratum rather than
combining strata in one LUND file.  Final campaign weights are produced with:

```tcsh
python3 radiative_full_campaign.py finalize \
  --campaign "$pilots/campaign.json"
```

The resulting `campaign_weights.tsv` contains one pooled cross section and
event weight per stratum while preserving the per-replica provenance needed
to connect each LUND file to downstream simulation.
