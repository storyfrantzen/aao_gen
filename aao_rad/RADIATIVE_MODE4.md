# Development radiative mode 4

Radiative mode 4 calibrates and then generates unweighted events for one
final-LUND analysis stratum per invocation. It is a development prototype:
its proposal and normalization are exact, but the current guard recipes and
tuning have not yet passed fresh-replica production validation.

The target coordinates are reconstructed from the final LUND electron and
proton after internal radiation and both incoming and outgoing external
energy loss:

```text
Q2_observed, xB_observed, minus_t_observed, phi_observed_deg
```

The configured observed `W_min` is applied. No lower `y` cut exists, and the
configured upper `y` cut is disabled unless `--apply-y-max` is requested
explicitly.

## Exact proposal

For stratum `i`, the proposal is

```text
q_i = alpha * q_core,i + (1 - alpha) * q_legacy
```

where `q_core,i` is uniform in the learned five-dimensional continuous box:

```text
r_u
r_ep
u_gamma
hadron_cosine_base
hadron_phi_base
```

The remaining legacy random variables, including radiative channel,
photon-angle, target-loss, and rotation variables, are unrestricted.
`hadron_phi_base` uses the recipe's periodic origin and wrapped interval.
Every reconstructed guard is anchored at `u_gamma = 1`, the soft-photon
endpoint. The lower `u_gamma` face still comes from the learned recipe and
padding. The anchored bounds, their newly computed volume, and the anchor
flag are frozen in the manifest. This prevents a harmless soft-endpoint
sample from being treated as a high-weight outside-tail event merely because
the finite survey did not reach exactly to one.

The legacy component has strictly positive probability. Consequently, a
guard that misses a physical region can reduce efficiency but cannot bias the
generated distribution. Every internally valid proposal is corrected by

```text
q_legacy / q_i =
    1 / ((1 - alpha) + alpha / V_i)   inside the core
    1 / (1 - alpha)                   outside the core
```

where `V_i` is the normalized five-dimensional core volume. The final
candidate is constructed and classified before acceptance-rejection. This is
required because outgoing external radiation can migrate an event between
final-LUND strata.

AAO retains its exact stochastic multiplicity correction when the supplied
`sigr_max` is exceeded:

```text
mcall = floor(sigr / sigr_max)
      + Bernoulli(frac(sigr / sigr_max))
```

Thus `mcall_max > 1` does not bias the distribution, but it creates duplicate
complete candidates and signals that a larger envelope may be operationally
preferable.

## Evidence-based face refinements

A finite training survey can stop just short of an important continuous-guard
face. Do not edit the learned recipe or generated AAO input by hand. Mode 4
can instead create and apply an optional evidence-hashed refinement artifact.
Refinement values are the exact **final post-padding guard bounds**. They are
not padded again, and they may expand support but may never contract it.

For the first `s04468` refinement:

```bash
python3 radiative_mode4.py create-refinement \
  --config ../../../configs/analysis/rgk/6.535.json \
  --recipes \
    rgk_mode4_v3_fresh_20260730/inputs/continuous_guard_recipes.json \
  --candidate padding_0p035 \
  --output \
    rgk_mode4_v3_fresh_20260730/inputs/s04468_refinement_iteration001.json \
  --stratum s04468 \
  --face r_u:upper:0.670 \
  --face r_ep:lower:0.275 \
  --rationale \
    "Independent complement calibration found two r_u-upper escapes and one r_ep-lower escape." \
  --evidence \
    rgk_mode4_v3_fresh_20260730/calibration/s04468_pooled.json \
  --evidence \
    rgk_mode4_v3_fresh_20260730/calibration/s04468_complement_5M/runs/s04468/s04468__g0000.calibration.csv
```

The creator validates the analysis and recipe hashes, the guard candidate,
axis names, native coordinate limits, endpoint anchor, and expansion
direction. It records hashes and byte counts for every evidence artifact and
previews the original and refined boxes plus their volume ratio.

Pass the resulting artifact to either preparation operation:

```bash
python3 radiative_mode4.py prepare-calibration \
  --config ../../../configs/analysis/rgk/6.535.json \
  --recipes \
    rgk_mode4_v3_fresh_20260730/inputs/continuous_guard_recipes.json \
  --refinements \
    rgk_mode4_v3_fresh_20260730/inputs/s04468_refinement_iteration001.json \
  --input aao_input.inp \
  --output \
    rgk_mode4_v3_fresh_20260730/calibration/s04468_refined_iteration001 \
  --candidate padding_0p035 \
  --core-fraction 0.90 \
  --inside-guard-trial-fraction 0.50 \
  --trials 1000000 \
  --heartbeat-interval 10000 \
  --replicas 1 \
  --bin-start 4468 \
  --bin-stop 4469 \
  --generator-revision `git rev-parse HEAD`
```

Preparation snapshots the complete refinement JSON. The manifest records its
SHA-256, lists the named refined strata, and freezes, for every run:

- the original learned-and-padded guard;
- the exact refined guard;
- every changed face and signed movement;
- the rationale and evidence hashes;
- the original and refined volumes and their ratio.

Run also verifies the frozen refinement snapshot before invoking AAO, and the
hash is copied into the normalization output and downstream summary.
Calibration campaigns can be pooled only when their refinement hashes and
final guard bounds agree. In particular, do not pool pre-refinement and
post-refinement `s04468` results.

## Calibrate the envelope first

Build from `aao_rad`:

```bash
make
make test
```

Do not guess `sigr_max`. A value large enough for the amplified legacy tail
can make the densely sampled core extraordinarily inefficient, while a value
near the core scale can make a rare tail proposal emit a large multiplicity.
The fixed-trial calibration operation samples two **disjoint regions**
without acceptance-rejection or LUND output:

```text
C       = the learned guard
C-bar   = the exact complement of that guard in the native unit hypercube
```

The inside component is uniform on `C`. The complement component is sampled
uniformly on `C-bar` by rejecting the guard itself; therefore every noncore
calibration trial is known to be outside the guard. This is a calibration
proposal only. Production remains the 90% learned-core plus 10%
unrestricted-legacy full-support mixture described above.

If the production core fraction is `alpha` and the guard volume is `V`, the
production proposal puts these probability masses in the two calibration
regions:

```text
P(C)     = alpha + (1 - alpha) V
P(C-bar) = (1 - alpha) (1 - V)
```

If calibration assigns a fraction `beta` of trials inside the guard, its
regional importance factors are:

```text
inside guard:      P(C)     / beta
guard complement:  P(C-bar) / (1 - beta)
```

These factors preserve both the intended production proposal and the
physical cross-section integral. `beta` affects precision and runtime only;
it does not alter the result. Campaigns with different values of `beta` can
therefore be pooled component by component.

Prepare a 100,000-trial first calibration for one representative stratum:

```bash
python3 radiative_mode4.py prepare-calibration \
  --config ../../../configs/analysis/rgk/6.535.json \
  --recipes \
    migration_rgk_balanced_continuous_guards_iteration001/continuous_guard_recipes.json \
  --input aao_input.inp \
  --output mode4_rgk_calibration_s04468 \
  --candidate padding_0p035 \
  --core-fraction 0.90 \
  --inside-guard-trial-fraction 0.50 \
  --trials 100000 \
  --heartbeat-interval 10000 \
  --replicas 1 \
  --bin-start 4468 \
  --bin-stop 4469 \
  --generator-revision `git rev-parse HEAD`
```

The input file must be a legacy AAO input without a survey trailer and must
match the analysis beam energy. The wrapper freezes `fmcall=0`; calibration
uses an unused positive placeholder envelope, while generation freezes the
explicit `--sigr-max`. The original file is not modified.

The default does not apply an analysis-level upper-`y` cut. Add
`--apply-y-max` only after that selection is intentionally frozen.

Run the fixed-trial calibration:

```bash
python3 radiative_mode4.py run \
  mode4_rgk_calibration_s04468/manifest.json \
  --flat-index 4468 \
  --replica-index 0 \
  --executable build/aao_rad \
  |& tee mode4_rgk_calibration_s04468_run.log
```

At startup the wrapper prints `mode4_live_heartbeat=...`. Inspect it safely
from another session:

```bash
tail -f /printed/path/aao_rad.mode4.heartbeat.csv
```

The heartbeat is flushed at the requested proposal interval and contains
proposal, internal-valid, final-candidate, target-candidate, event, and
inside/noncore component counts.

Finalize the calibration:

```bash
python3 radiative_mode4.py finalize \
  mode4_rgk_calibration_s04468/manifest.json \
  --envelope-safety-factor 1.20 \
  --maximum-duplicate-fraction 0.05 \
  --minimum-component-targets 20
```

The resulting `envelope_calibration.json` reports:

- inside-guard and guard-complement target rates, cross sections,
  uncertainties, and ESS;
- the exact production probability mass assigned to each disjoint region;
- corrected-integrand quantiles and observed maxima;
- expected event yield, emitting-proposal rate, duplicate fraction, and
  maximum observed `mcall` ratio for each candidate envelope;
- a recommendation only when both disjoint regions have enough target
  support.

If the report says `insufficient_guard_complement_target_support`, add a
complement-heavy calibration rather than discarding the completed work:

```bash
python3 radiative_mode4.py prepare-calibration \
  --config ../../../configs/analysis/rgk/6.535.json \
  --recipes \
    migration_rgk_balanced_continuous_guards_iteration001/continuous_guard_recipes.json \
  --input aao_input.inp \
  --output mode4_rgk_calibration_s04468_complement_heavy \
  --candidate padding_0p035 \
  --core-fraction 0.90 \
  --inside-guard-trial-fraction 0.20 \
  --trials 1000000 \
  --heartbeat-interval 10000 \
  --replicas 1 \
  --bin-start 4468 \
  --bin-stop 4469 \
  --seed-base 682001 \
  --generator-revision `git rev-parse HEAD`
```

After running that manifest, pool both campaigns in one finalization:

```bash
python3 radiative_mode4.py finalize \
  mode4_rgk_calibration_s04468/manifest.json \
  mode4_rgk_calibration_s04468_complement_heavy/manifest.json \
  --envelope-safety-factor 1.20 \
  --maximum-duplicate-fraction 0.05 \
  --minimum-component-targets 20 \
  --output mode4_rgk_calibration_s04468_pooled.json
```

Pooling requires identical generator revision, configuration and recipe
hashes, guard candidate, production core fraction, selections, and per-stratum
guard bounds. Trial allocation, trial count, replica count, and seed may
differ. Duplicate manifests and repeated stratum seeds are rejected.

## Generate after calibration

Use the recommended `sigr_max` from the frozen calibration report:

```bash
python3 radiative_mode4.py prepare \
  --config ../../../configs/analysis/rgk/6.535.json \
  --recipes \
    migration_rgk_balanced_continuous_guards_iteration001/continuous_guard_recipes.json \
  --input aao_input.inp \
  --output mode4_rgk_pilot_s04468 \
  --candidate padding_0p035 \
  --core-fraction 0.90 \
  --sigr-max RECOMMENDED_VALUE \
  --events-per-stratum 5000 \
  --heartbeat-interval 100000 \
  --replicas 1 \
  --bin-start 4468 \
  --bin-stop 4469 \
  --generator-revision `git rev-parse HEAD`
```

Run and finalize generation:

```bash
python3 radiative_mode4.py run \
  mode4_rgk_pilot_s04468/manifest.json \
  --flat-index 4468 \
  --replica-index 0 \
  --executable build/aao_rad

python3 radiative_mode4.py finalize \
  mode4_rgk_pilot_s04468/manifest.json
```

`finalize` pools replica cross sections using their proposal counts:

```text
combined_sigma_i =
    sum_r(ntries_ir * sigma_ir) / sum_r(ntries_ir)

event_weight_i = combined_sigma_i / sum_r(Nevents_ir)
```

Those are the weights for combining strata downstream. Within one stratum,
the LUND events are unweighted and distributed according to the physical
radiative cross section conditioned on that final-LUND stratum.

## Artifacts and type-2 compatibility

`prepare` snapshots the exact analysis configuration, continuous recipes, and
legacy input. Its immutable manifest records:

- canonical `sNNNNN` stratum and `gNNNN` generation identifiers;
- requested bin bounds and indices;
- generator, config, and recipe revisions or hashes;
- reconstructed padded guard and normalized volume;
- optional evidence-hashed guard refinements, with original and refined
  bounds;
- core and unrestricted-tail fractions;
- explicit seed, event count, and `sigr_max`.

`run` writes products under a stem such as:

```text
runs/s04468/s04468__g0000
```

The `.norm` and run JSON record the stratum cross section, event weight,
proposal count, core/noncore trial and event counts, maximum multiplicity,
and proposal efficiency. `aao_rad.mode4.csv` records accepted-event diagnostics
and is independently checked against the manifest proposal density and
final-coordinate bounds.

Every run also preserves `*.heartbeat.csv`. Calibration runs preserve
`*.calibration.csv`, emit no LUND events, and are summarized in
`envelope_calibration.json` plus a compact TSV.

The LUND structure is unchanged: one historical event header followed by the
same four particle records. No weights or metadata are inserted into LUND.
The stratum ID remains in the filename and manifest, making the products
compatible with OSG type-2 LUND submission when each chunk contains only one
stratum.

## Development validation

Start with a small representative set rather than all analysis strata:

- high-cross-section, well-covered strata;
- low- and high-kinematic strata;
- low-coverage cases such as `s04468` and `s09012`.

Calibrate and then compare `padding_0p03`, `padding_0p035`, and
`padding_0p04`, initially at a 90/10 generation mixture. Inspect:

- emitted events per proposal;
- `mcall_max` and duplicate-event frequency;
- core-versus-legacy emitted fractions;
- stratum-integrated cross-section stability across replicas;
- consistency with weighted survey estimates;
- exact final-LUND coordinate membership.

Once the algorithm and tuning are frozen, refit recipes using all development
replicas and validate the frozen result on fresh survey replicas. A global
radiative sample remains necessary for closure.

## Tests

Run the focused mode-4 tests:

```bash
python3 -m unittest -v test_radiative_mode4.py
```

Run the complete radiative workflow suite:

```bash
make test
```
