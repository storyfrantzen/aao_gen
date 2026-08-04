# Active radiative-stratum workflow

`radiative_active_strata.py` separates three questions that must not be
answered from the same Monte Carlo count:

1. Does the nominal analysis hypercube intersect physical phase space?
2. Can the stratum affect the measurement directly or through migration?
3. Is its mode-4 calibration and generation validation ready?

A zero survey or calibration count never proves structural emptiness. A
stratum may be marked `structurally_empty` only by explicit, cited external
evidence. Similarly, `closure_only` requires complete evidence that the bin
is not measured and that its model cross section, feed-in, and global closure
impact all lie below declared analysis thresholds.

## Classification states

| Status | Meaning |
|---|---|
| `active_ready` | Relevant, calibrated, and pilot-validated. |
| `active_needs_pilot_validation` | Relevant with a usable envelope, but independent generation pilots are incomplete. |
| `active_needs_optimization` | Relevant but its calibration/envelope is not ready. |
| `closure_only` | Nonempty but demonstrably negligible for the selected analysis thresholds; retain it in global mode-3 closure. |
| `structurally_empty` | Explicit evidence proves the physical intersection is empty. |
| `needs_relevance_assessment` | Evidence is missing or insufficient; no exclusion decision is allowed. |

The first three statuses form `active_flat_indices.txt`. Only
`active_ready` appears in `production_ready_flat_indices.txt`.

## Build reconstructed-data occupancy evidence

Build the full analysis-catalog occupancy table directly from the selected
data artifact and its exclusivity mask:

```bash
python3 radiative_active_strata.py build-data-evidence \
  --config analysis_ymax0p95.json \
  --data-events data_events.npz \
  --selection-mask data_exclusivity.npy \
  --output data_occupancy_iteration000
```

The command reapplies the configured `Q2_min`, `W_min`, and optional `y_max`
selection in reconstructed coordinates before assigning the half-open
analysis bins. It verifies array lengths, requires a boolean selection mask,
checks selected `(run,event)` keys for accidental duplication, and records
SHA-256 digests for every input. Duplicate keys fail by default; use
`--allow-duplicate-event-keys` only after demonstrating that the overlap is
intentional.

The immutable output directory contains:

- `data_occupancy.json`: cut flow, coordinate definitions, input provenance,
  catalog summary, and all per-stratum counts;
- `data_occupancy.tsv`: a compact 12,960-row review table for the RGK binning;
- `relevance_evidence.json`: a directly compatible input to `classify`;
- `occupied_flat_indices.txt` plus lists requiring at least 5, 10, or 50 data
  events.

Positive data occupancy makes a stratum relevant to the analysis, but the
builder deliberately leaves `physical_status` and `analysis_included`
unknown. A reconstructed event does not by itself prove generator-level
physical support, and a zero count never proves structural emptiness or
negligibility. Model, feed-in, and global-closure evidence must be added
separately before an unoccupied stratum can become `closure_only`.

## Add pooled survey-model evidence

After freezing a migration training manifest and its independent holdout
validation, augment the data relevance without overwriting it:

```bash
python3 radiative_active_strata.py augment-survey-evidence \
  --config analysis_ymax0p95.json \
  --base-relevance data_occupancy_iteration000/relevance_evidence.json \
  --migration-manifest migration_iteration000/migration_manifest.json \
  --migration-validation \
    migration_iteration000_validation/migration_validation.json \
  --output survey_relevance_iteration000
```

The command requires a full-catalog base relevance artifact. It verifies all
configuration and manifest hashes and requires the migration `Q2`, `W`, and
optional `y` policy to match the analysis configuration exactly. Training and
holdout sufficient statistics are pooled as one fixed-trial estimator:

```text
pooled cross section = (training contribution sum + holdout contribution sum)
                       / (training proposals + holdout proposals)
```

The sum-of-squares and proposal counts are pooled at the same time, so the
reported SEM and importance ESS are not averages of the two replica metrics.
Each stratum's `model_cross_section_fraction` uses the pooled inside-analysis
cross section as its denominator, and the builder verifies that all 12,960
fractions close to one.

The immutable output contains:

- `survey_model_evidence.json`: full provenance, pooling closure, support and
  parent-coverage counts, and per-stratum training/holdout/pooled metrics;
- `survey_model_evidence.tsv`: the compact full-catalog review table;
- `relevance_evidence.json`: the data evidence augmented with model fractions;
- lists for nonzero survey support, independent support, and failed frozen
  parent coverage.

No survey contribution is interpreted as zero physical cross section.
Likewise, failed hard-parent coverage means that a guard needs work; it never
removes an otherwise relevant stratum. This artifact does not populate
`maximum_feed_in_fraction`, because detector feed-in requires a separately
audited GEMC response.

## Build a cumulative calibration queue

A per-stratum model threshold can discard a collectively important diffuse
tail. Select every data-occupied stratum and then add zero-data strata in
descending pooled-model order until the remaining global model fraction is
below an explicit residual budget:

```bash
python3 radiative_active_strata.py build-cumulative-queue \
  --config analysis_ymax0p95.json \
  --relevance survey_relevance_iteration000/relevance_evidence.json \
  --minimum-data-events 1 \
  --maximum-global-model-residual-fraction 0.001 \
  --output cumulative_queue_iteration000
```

Ties in model fraction are resolved by ascending flat index. Data strata are
ordered first by descending event count, then model fraction, so the output is
both deterministic and useful as a calibration priority list. The builder
requires complete data counts, complete pooled model fractions that close to
one, and the survey support metadata written by `augment-survey-evidence`.

Selected strata are routed without changing their relevance decision:

| Work category | Meaning |
|---|---|
| `supported_calibration` | Independent survey support and passed frozen-parent coverage; start calibration directly. |
| `guard_refinement` | Some survey support exists, but it is one-sided or parent coverage failed; improve the guard evidence first. |
| `targeted_discovery` | No pooled survey contribution was observed; run targeted discovery rather than declaring the stratum empty. |

The immutable output contains:

- `cumulative_stratum_queue.json` with configuration/evidence hashes, exact
  policy, achieved residual, work-category counts, and all catalog records;
- `cumulative_stratum_queue.tsv`, ordered with selected work first;
- `selected_flat_indices.txt`, directly usable as a mode-4 campaign selector;
- separate data-occupied and model-required-zero-data lists;
- separate supported-calibration, guard-refinement, and targeted-discovery
  work lists;
- `omitted_flat_indices.txt`, retained for global mode-3 and later feed-in
  closure checks.

Omitted means only that the current pooled model contribution lies within the
declared global residual. It does not establish structural emptiness. GEMC
feed-in evidence can reactivate an omitted stratum in a later queue revision.

## Create an evidence template

Create the template only after identifying the analysis artifacts that will
support its claims. Every cited artifact is stored with an absolute path and
SHA-256 digest:

```bash
python3 radiative_active_strata.py create-template \
  --config analysis_config.json \
  --output relevance_evidence.json \
  --flat-index 10963 \
  --evidence-source data=data_occupancy.json \
  --evidence-source migrations=mode3_migration_relevance.json \
  --evidence-source closure=mode3_closure.json
```

Without `--flat-index`, the template contains the complete analysis catalog.
Each stratum starts with unknown fields:

```json
{
  "physical_status": "unknown",
  "analysis_included": null,
  "data_events": null,
  "model_cross_section_fraction": null,
  "maximum_feed_in_fraction": null,
  "global_closure_impact_fraction": null,
  "force_active": false,
  "rationale": null,
  "source_ids": []
}
```

Allowed physical states are `unknown`, `nonempty`, `partially_accessible`,
and `structurally_empty`. Any nondefault claim requires a nonempty rationale
and at least one valid `source_id`. Fractions are defined relative to the
full analysis prediction, not merely the subset of calibrated mode-4 bins.

Positive fixed-trial calibration evidence is sufficient to establish that a
stratum is nonempty. A zero calibration is deliberately not sufficient to
establish that it is empty.

## Classify

Thresholds are mandatory because the generator must not invent an analysis
definition of “negligible”:

```bash
python3 radiative_active_strata.py classify \
  --config analysis_config.json \
  --relevance relevance_evidence.json \
  --calibration canonical_envelopes.json \
  --pilot-validation canonical_pilot_validation.json \
  --flat-index 10963 \
  --minimum-data-events 1 \
  --maximum-model-fraction 0.001 \
  --maximum-feed-in-fraction 0.001 \
  --maximum-closure-fraction 0.001 \
  --output active_strata_iteration000
```

Multiple disjoint canonical calibration and pilot-validation reports may be
supplied by repeating their options. Overlapping reports are rejected so a
stale calibration cannot silently override the report selected for a
stratum. A pilot-validation report is accepted only when its referenced
calibration hash is among the supplied canonical calibrations.

If an audited full-model integrated cross section is available, pass
`--reference-cross-section-microbarn VALUE`. A missing per-stratum model
fraction is then derived from its calibration cross section and this
denominator. Data occupancy, feed-in, and global closure impact remain
independent required evidence before a stratum can become `closure_only`.

## Outputs

Classification creates an immutable directory containing:

- `active_stratum_mask.json`: complete evidence, thresholds, provenance, and
  decisions;
- `active_stratum_mask.tsv`: compact review table;
- `active_flat_indices.txt` and `production_ready_flat_indices.txt`;
- one flat-index list for every classification status.

These lists are campaign selectors, not physics weights. Each generated
stratum retains its own calibrated cross section and event weight.

Both `radiative_mode4.py prepare` and `prepare-calibration` accept a selected
list directly:

```bash
python3 radiative_mode4.py prepare \
  --flat-index-file active_strata_iteration000/production_ready_flat_indices.txt \
  ...
```

The selection file must contain one integer per line; blank lines and lines
beginning with `#` are ignored. It is mutually exclusive with `--flat-index`,
`--bin-start`, and `--bin-stop`. Preparation copies the file into the campaign,
records its SHA-256 in the manifest, and execution rejects a changed snapshot.

## Required closure check

`closure_only` means “do not optimize with bin-conditional mode 4,” not
“delete from physics.” Closure-only strata remain represented by the global
mode-3 campaign. Before production, compare reconstructed yields and the
response matrix with and without their explicit mode-4 columns. The observed
change must remain below the same closure and migration thresholds stored in
the mask.
