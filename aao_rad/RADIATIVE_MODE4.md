# Development radiative mode 4

For deciding which nominal analysis strata should receive bin-conditional
optimization, see [ACTIVE_STRATUM_WORKFLOW.md](ACTIVE_STRATUM_WORKFLOW.md).
That workflow separates physical support, measurement/migration relevance,
and generator readiness; zero Monte Carlo support never silently removes a
stratum.

To quantify how an already calibrated loose campaign survives a stricter
final-electron momentum cut, see
[RADIATIVE_MOMENTUM_AUDIT.md](RADIATIVE_MOMENTUM_AUDIT.md).

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
inside/noncore component counts. Production heartbeats also record the
number of distinct emitting candidates and the exact number of additional
events created by `mcall > 1`; their ratio to all events is the observed
duplicate-event fraction.

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
- a strict recommendation only when both disjoint regions have enough target
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

If generator revisions differ only because of audited diagnostic or
finalizer changes, calibration pooling can be explicitly authorized with:

```bash
  --allow-calibration-revision-mismatch \
  --revision-compatibility-rationale \
    "Exact reason the calibration physics and proposal are unchanged"
```

This override is calibration-only, is rejected without a nonempty rationale,
and records every source revision plus the rationale in the report. It must
not be used across changes to the integrand, phase-space mapping, proposal,
selection, or guard construction.

### Provisional zero-complement stopping rule

A well-tuned guard may produce no complement targets even after a large,
independent complement-directed campaign. The default finalizer remains
strict and will not recommend an envelope in that case. An explicit opt-in
policy can instead certify a small production pilot:

```bash
python3 radiative_mode4.py finalize \
  refined_equal_allocation/manifest.json \
  refined_complement_heavy/manifest.json \
  --envelope-safety-factor 1.20 \
  --maximum-duplicate-fraction 0.05 \
  --minimum-component-targets 20 \
  --minimum-provisional-inside-targets 1000 \
  --additional-inside-pilot-run \
    previous_pilot/runs/s04468/s04468__g0000.json \
  --allow-zero-complement \
  --zero-complement-confidence 0.95 \
  --maximum-zero-complement-target-rate 1e-6 \
  --output refined_provisional_envelope.json
```

For zero observed complement targets in `N` independent complement trials,
the report evaluates the exact one-sided binomial upper occurrence rate:

```text
p_upper = 1 - (1 - confidence)^(1/N)
```

The provisional policy is allowed only when:

- the inside guard meets `--minimum-component-targets`;
- at least `--minimum-provisional-inside-targets` inside-guard target
  candidates have been observed (default 1000);
- exactly zero complement targets were observed;
- `p_upper` does not exceed the configured maximum target rate;
- the safety-scaled maximum corrected integrand observed inside the guard
  meets the duplicate-fraction requirement.

Each optional `--additional-inside-pilot-run` must be a validated mode-4
generation run containing no noncore events and matching the finalized
stratum's frozen bounds and guard. Its event-file maximum can only raise the
envelope floor. The report records hashes of both artifacts, and the pilot
events are deliberately excluded from fixed-trial cross-section, uncertainty,
ESS, occurrence-rate, and target-count calculations.

If even one complement target is observed, the zero-target exception is
disabled and the ordinary component-support requirement remains in force.
The resulting status is `provisional_zero_complement`, not `recommended`.
The JSON and TSV record:

- the confidence level, complement trial count, exact upper occurrence rate,
  and configured threshold;
- the finalizer's repository revision and exact source-file hash;
- the recommendation basis and provisional flag;
- the inside-target threshold, observed count, empirical next-target rank
  resolution, observed maximum, and safety-scaled maximum;
- any inside-only pilot evidence, its artifact hashes, and the additional
  observed maximum used as an envelope floor;
- `pilot_readiness=ready_provisional_zero_complement`;
- a warning that the occurrence-rate bound does **not** bound an unseen
  event's cross section or corrected-integrand magnitude.

The empirical next-target rank resolution, `1/(n_inside+1)`, is an audit
indicator rather than a confidence bound on the integrand tail. The envelope's
predicted yield and duplicate metrics remain conditional on the observed
calibration sample. Full proposal support and stochastic multiplicity preserve
correctness if a later event exceeds the envelope, but such an event can
create duplicates and reduce effective precision. Therefore this policy
authorizes only a small monitored pilot before broader production.

## Generate after calibration

For a one-stratum or deliberately conservative smoke test, a shared scalar
envelope remains supported:

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

For a heterogeneous multi-stratum campaign, use the finalized calibration
report directly instead. Each selected stratum must have
`pilot_readiness=ready` or `ready_provisional_zero_complement` and a positive
`recommended_envelope.sigr_max`:

```bash
python3 radiative_mode4.py prepare \
  --config ../../../configs/analysis/rgk/6.535.json \
  --recipes continuous_guard_recipes.json \
  --refinements guard_refinements.json \
  --input aao_input.inp \
  --output mode4_rgk_sparse_pilot \
  --candidate padding_0p035 \
  --core-fraction 0.90 \
  --envelope-report multi_stratum_envelope_calibration.json \
  --events-per-stratum 200 \
  --heartbeat-interval 100000 \
  --replicas 2 \
  --flat-index 4468 \
  --flat-index 9012 \
  --seed-base 997101 \
  --generator-revision `git rev-parse HEAD`
```

`--flat-index` is repeatable and selects disjoint strata without preparing
the bins between them. It is mutually exclusive with `--bin-start` and
`--bin-stop`. A contiguous range may still be used with an envelope report.
An audited active-mask list can instead be supplied with
`--flat-index-file`; it is snapshotted and hashed by the prepared manifest.

The report must match the selected analysis-config, recipe, refinement,
guard-candidate, core-fraction, analysis-selection, bounds, and reconstructed
guard hashes and values. A selected missing or unready stratum is rejected.
If a wrapper-only revision separates the calibration and generation commits,
an explicit audited exception is required:

```bash
  --allow-envelope-revision-mismatch \
  --envelope-revision-compatibility-rationale \
    "Wrapper-only change; Fortran physics and proposal are unchanged."
```

The complete calibration report is copied to
`envelope_calibration.json` inside the prepared campaign and hashed by the
manifest. Every run record freezes its own `sigr_max`, readiness state, and
recommendation basis. Execution and pilot validation reject a changed
snapshot or a run-level envelope inconsistent with the manifest. Legacy v3
shared-envelope manifests remain readable.

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

## Validate independent pilots

Before expanding a calibrated stratum to a broader campaign, pool at least
two independent generation pilots with `validate-pilots`:

```bash
python3 radiative_mode4.py validate-pilots \
  --calibration s04468_envelope_calibration.json \
  --run pilot_1/runs/s04468/s04468__g0000.json \
  --run pilot_2/runs/s04468/s04468__g0000.json \
  --minimum-runs 2 \
  --minimum-events 400 \
  --maximum-duplicate-fraction 0.05 \
  --maximum-guard-complement-fraction 0.02 \
  --maximum-relative-cross-section-difference 0.10 \
  --maximum-cross-section-z-score 3 \
  --confidence 0.95 \
  --output s04468_pilot_validation.json
```

The command verifies every run and its hashed source manifest and event CSV,
rejects repeated seeds, and requires matching analysis, recipe, refinement,
guard, selection, audited generator revision, and stratum metadata. A
revision mismatch is rejected unless
`--allow-pilot-revision-mismatch` is accompanied by a nonempty
`--revision-compatibility-rationale`; the override is preserved in the
report. It writes JSON plus a compact TSV. For each stratum it reports:

- proposal-count-weighted pilot cross section and run-to-run uncertainty;
- closure against the fixed-trial calibration cross section;
- exact emitting-candidate and multiplicity-generated duplicate counts;
- duplicate overhead and a one-sided Wilson diagnostic;
- event and emitting-candidate counts in the four combinations of
  guard-focused/legacy proposal component and inside/outside geometric guard;
- maximum corrected integrands by proposal component and geometric region;
- every crossed guard face and its largest observed excursion;
- `hold` versus `increase_and_revalidate` envelope guidance;
- `hold` versus `review_complement_geometry` guard guidance;
- a final `ready_for_multi_stratum_pilot` recommendation.

Proposal component and geometric membership are deliberately kept separate.
The unrestricted legacy component can land either inside or outside the
guard. A guard-complement event with stochastic multiplicity does not by
itself fail validation: measured duplicate overhead, complement frequency,
and cross-section closure determine the decision.

Pilots generated with an envelope no larger than the currently recommended
envelope are conservative stress tests for that target and count toward its
support. A rare integrand above `sigr_max` is recorded, but does not force an
envelope increase when stochastic multiplicity keeps the pooled duplicate
overhead below threshold. Wilson bounds are useful finite-sample diagnostics,
not formal guarantees for nonidentically distributed candidates.

The default readiness thresholds require two independent runs and 400 total
qualifying events. These defaults are intended for a development pilot, not
final production certification.

## Artifacts and type-2 compatibility

`prepare` snapshots the exact analysis configuration, continuous recipes, and
legacy input. When used, it also snapshots the finalized envelope-calibration
report. Its immutable v4 manifest records:

- canonical `sNNNNN` stratum and `gNNNN` generation identifiers;
- requested bin bounds and indices;
- generator, config, and recipe revisions or hashes;
- reconstructed padded guard and normalized volume;
- optional evidence-hashed guard refinements, with original and refined
  bounds;
- core and unrestricted-tail fractions;
- explicit seed, event count, and per-run `sigr_max`;
- contiguous-range or sparse-flat-index selection;
- envelope mode, report SHA-256, selected readiness states, and any explicit
  revision-compatibility audit.

`run` writes products under a stem such as:

```text
runs/s04468/s04468__g0000
```

The `.norm` and run JSON record the resolved per-stratum envelope, stratum
cross section, event weight,
proposal count, core/noncore trial and event counts, maximum multiplicity,
distinct emitting candidates, exact multiplicity-generated duplicate count
and fraction, and proposal efficiency. `aao_rad.mode4.csv` records
accepted-event diagnostics and is independently checked against the manifest
proposal density and final-coordinate bounds.

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
