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
