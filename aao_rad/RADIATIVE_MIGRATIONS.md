# Milestone 2b hard-parent migration diagnostic

`radiative_migrations.py` tests a physically motivated middle ground between
the exact Born analysis strata and an independent high-dimensional
random-number guard for every radiative stratum. It measures the
cross-section-weighted mapping

```text
hard/Born-like parent + radiative channel -> final-LUND analysis stratum
```

This is a diagnostic workflow. It does not change AAO sampling, implement mode
4, perform acceptance-rejection, or write LUND events.

## Coordinate definitions

The final target stratum uses the existing `final_lund_analysis` definition:

```text
q2_observed, xb_observed, minus_t_observed, phi_observed_deg
```

The hard parent uses:

```text
q2_hard, xb_hard, minus_t_hard, phi_cm_deg
```

The hard coordinates are assigned with the same `(Q2, xB, -t, phi)` edges as
the final analysis strata. Hard `Q2`, `xB`, and `-t` values below or above
those edges are retained in explicit underflow and overflow parents. Hard
`phi` is periodic. `intreg` remains a separate radiative-channel component.

The final-stratum assignment applies the configured observed `Q2` and `W`
selection. It applies no analysis-level `y` cut unless `--apply-y-max` is
explicitly passed while learning. Validation inherits the frozen policy.

## Learn iteration zero

The five existing RGK pilot replicas are sufficient for the first diagnostic.
From `aao_rad`, use replicas 0, 1, and 2 for training:

```bash
python3 radiative_migrations.py learn \
  --config ../../../configs/analysis/rgk/6.535.json \
  --survey \
    survey_rgk_replica000 \
    survey_rgk_replica001 \
    survey_rgk_replica002 \
  --target-parent-fraction 0.995 \
  --parent-dilation 0 \
  --generator-revision 54f5ba8d59b86b8cabb2229e5c2cf2be5de1ff00 \
  --output migration_rgk_iteration000
```

For each final analysis stratum, occupied hard-parent/channel components are
ranked by their summed fixed-trial cross-section contribution. The smallest
set reaching `--target-parent-fraction` becomes the seed footprint.
`--parent-dilation` adds repeated axis-neighbor layers to that footprint;
hard phi wraps periodically and dilation never crosses radiative channels.

The default dilation is zero so that the diagnostic first measures the
compact Born-parent footprint itself. Every output also reports what one
additional dilation layer would recover.

## Validate on independent replicas

Use replicas 3 and 4 only after the parent footprint is frozen:

```bash
python3 radiative_migrations.py validate \
  --manifest migration_rgk_iteration000/migration_manifest.json \
  --survey \
    survey_rgk_replica003 \
    survey_rgk_replica004 \
  --minimum-parent-coverage 0.98 \
  --output migration_rgk_iteration000_validation
```

The command exits with status 2 when any assessed stratum fails the requested
held-out coverage. That is a diagnostic result, not a corrupt artifact.

Coverage is cross-section weighted:

```text
held-out target contribution from selected hard parents
--------------------------------------------------------
         all held-out contribution in the target
```

The reported cross-section purity is a different quantity:

```text
target-stratum contribution from selected hard parents
-------------------------------------------------------
 global final-valid contribution from those same parents
```

Purity estimates how selectively a parent footprint identifies one target
stratum under the surveyed legacy distribution. It is not yet a mode-4
acceptance efficiency because no hard-parent proposal has been implemented.

## Summarize and plot

Print the compact weighted summary:

```bash
python3 radiative_migrations.py summarize \
  --manifest migration_rgk_iteration000/migration_manifest.json \
  --validation \
    migration_rgk_iteration000_validation/migration_validation.json
```

Render the campaign diagnostics:

```bash
python3 radiative_migrations.py plot \
  --manifest migration_rgk_iteration000/migration_manifest.json \
  --validation \
    migration_rgk_iteration000_validation/migration_validation.json \
  --top-strata 6 \
  --output migration_rgk_iteration000_plots
```

The PDF contains cross-section-weighted coverage distributions, hard-parent
relationship fractions, aggregate hard-minus-observed bin offsets, and hard
`Q2` versus hard `xB` footprints for the highest-cross-section validation
strata.

## Artifacts

Learning writes:

- `migration_manifest.json`;
- `training_parent_coverage.csv`;
- `training_migrations.csv`.

Validation writes:

- `migration_validation.json`;
- `validation_parent_coverage.csv`;
- `validation_migrations.csv`.

Plotting writes:

- `migration_diagnostics.pdf`;
- `plot_summary.json`.

Every artifact receives a SHA-256 sidecar. Output directories are immutable
and must not already exist. Validation and plotting refuse to combine a
validation artifact with a different manifest hash.

The migration CSVs preserve, for every occupied
target/parent/channel component:

- final and hard bin indices;
- underflow/overflow labels;
- hard-minus-observed bin offsets;
- cross-section sum, squared sum, maximum, and contributing rows;
- fraction of the target-stratum cross section;
- seed, selected-parent, and one-more-dilation membership.

These sufficient statistics support later pooling and visualization without
rereading the full survey CSVs.

## Interpretation

A useful parent representation must balance:

- **completeness:** nearly all held-out target cross section is covered;
- **compactness:** relatively few parent components are selected;
- **purity:** the selected parents preferentially feed the target;
- **stability:** parent footprints persist across independent replicas.

The parent footprint is not a hard physics cut. A future mode-4 proposal must
retain a nonzero full-support tail and apply the exact proposal-density
correction before unweighting.
