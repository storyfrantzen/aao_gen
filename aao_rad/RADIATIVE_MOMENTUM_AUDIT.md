# Radiative final-electron momentum audit

`radiative_momentum_audit.py` measures how much of each loose radiative
analysis stratum is expected to survive a stricter final-electron momentum
selection.  It is intended for auditing an already learned mode-4 campaign
before deciding whether any strata need additional production optimization.
It does **not** change the guard geometry, rerun calibration, or introduce
event-by-event weights.

## Quantity being estimated

For observed stratum `i`, each fixed survey proposal has an importance-corrected
cross-section contribution `c_k`.  Missing CSV trials are zero contributions
and remain represented by the fixed proposal denominator `N`.

The loose and nominal estimates are

```text
sigma_loose_i   = sum(k in observed stratum i) c_k / N
sigma_nominal_i = sum(k in observed stratum i and p_e >= p_min) c_k / N

f_nominal_i = sigma_nominal_i / sigma_loose_i
```

Here `p_e` is calculated directly from the final LUND electron components:

```text
p_e = sqrt(final_e_px^2 + final_e_py^2 + final_e_pz^2)
```

Therefore `f_nominal_i` is the expected retained fraction of an unweighted
loose-stratum sample.  Its inverse is the expected number of loose generated
events needed per retained nominal event.  It is a production-efficiency
diagnostic, not a variable event weight.

The reported fraction SEM uses the delta method with the exact covariance of
the disjoint pass and fail contributions.  Raw row fractions are written only
as explicitly labelled diagnostics because balanced survey rows do not all
carry equal cross-section contribution.

## Recommended RGK command

Use the same frozen `y_max=0.95` analysis configuration that defined the
current survey/migration campaign.  The momentum threshold can be supplied
explicitly even if that frozen configuration predates `electron_p_min`:

```tcsh
set repo = /w/hallb-scshelf2102/clas12/storyf/SF_analysis_software_v2.0
set aao_rad = "$repo/external/aao_gen/aao_rad"
set fresh = "$aao_rad/rgk_mode4_v3_fresh_20260730"
set census = "$fresh/census/rgk_iteration000"
set queue = "$census/outputs/cumulative_queue_iteration000"

cd "$aao_rad"

# Set this to the frozen census config whose phase_space.y_max is 0.95.
set census_config = "$census/inputs/analysis_ymax0p95.json"
set momentum_audit = "$census/outputs/electron_pmin1_audit_iteration000"

python3 radiative_momentum_audit.py audit \
  --config "$census_config" \
  --training-survey "$fresh/surveys/replica005_10M" \
  --validation-survey "$fresh/surveys/replica006_10M" \
  --electron-p-min 1.0 \
  --apply-y-max \
  --flat-index-file "$queue/selected_flat_indices.txt" \
  --low-nominal-fraction 0.5 \
  --output "$momentum_audit"
```

Before starting, verify the intended denominator policy:

```tcsh
python3 -c 'import json,sys; p=json.load(open(sys.argv[1]))["phase_space"]; print("y_max =",p.get("y_max")); print("electron_p_min =",p.get("electron_p_min","supplied on command line"))' "$census_config"
```

If the surveys were not assigned training and validation roles, pool them with
one `--survey` option instead.  Do not combine `--survey` with the split
options.

## Artifacts

The output directory is immutable and contains:

- `electron_momentum_audit.json`: full provenance, estimator definition,
  selected-campaign summary, per-stratum sufficient statistics, and optional
  training/validation comparison;
- `electron_momentum_audit.tsv`: one row for every analysis stratum, with a
  `selected` flag for the 3,637-stratum campaign queue;
- `low_nominal_fraction_flat_indices.txt`: selected strata with observed
  `0 < f_nominal < 0.5` (or the requested planning threshold);
- `no_nominal_contribution_observed_flat_indices.txt`: selected strata with
  loose survey support but no observed nominal contribution;
- `unresolved_no_survey_contribution_flat_indices.txt`: selected strata for
  which the surveys cannot estimate a fraction;
- `selected_flat_indices.txt`: a normalized snapshot of the audited selection.

An observed fraction of zero or one is labelled as a finite-survey boundary;
it is not treated as proof that the physical fraction is exactly zero or one.

## Compact review commands

```tcsh
set audit_json = "$momentum_audit/electron_momentum_audit.json"
set audit_tsv = "$momentum_audit/electron_momentum_audit.tsv"

python3 -c 'import json,sys; r=json.load(open(sys.argv[1])); print(json.dumps(r["threshold"],indent=2)); print(json.dumps(r["summary"],indent=2))' "$audit_json"

awk -F '\t' 'NR==1 || ($3=="True" && $21!="" && $21+0<0.5)' \
  "$audit_tsv" | head -n 30

awk -F '\t' 'NR==1 || ($3=="True" && ($8=="no_survey_contribution" || $8=="no_nominal_contribution_observed"))' \
  "$audit_tsv" | head -n 30
```

The important planning columns are:

- `nominal_cross_section_fraction`: estimated survival probability for a
  loose mode-4 event from that stratum;
- `nominal_fraction_sem_delta_method`: survey uncertainty on that fraction;
- `loose_events_per_expected_nominal_event`: expected production overhead;
- `training_validation_fraction_difference_z_score`: an independent stability
  diagnostic when a split was supplied;
- denominator and nominal ESS: whether the estimate is dominated by a few
  importance-weighted survey trials.
