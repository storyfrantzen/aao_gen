# Development radiative mode 4

Radiative mode 4 generates unweighted events for one final-LUND analysis
stratum per invocation. It is a development prototype: its proposal and
normalization are exact, but the current guard recipes and tuning have not
yet passed fresh-replica production validation.

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

## Prepare a pilot

Build from `aao_rad`:

```bash
make
```

Prepare one representative stratum with the current development recipe:

```bash
python3 radiative_mode4.py prepare \
  --config ../../../configs/analysis/rgk/6.535.json \
  --recipes \
    migration_rgk_balanced_continuous_guards_iteration001/continuous_guard_recipes.json \
  --input aao_input.inp \
  --output mode4_rgk_pilot_s04468 \
  --candidate padding_0p035 \
  --core-fraction 0.90 \
  --sigr-max 0.05 \
  --events-per-stratum 5000 \
  --replicas 1 \
  --bin-start 4468 \
  --bin-stop 4469 \
  --generator-revision `git rev-parse HEAD`
```

The input file must be a legacy AAO input without a survey trailer and must
match the analysis beam energy. The wrapper freezes `fmcall=0` and the
explicit `--sigr-max`; the original file is not modified.

The default does not apply an analysis-level upper-`y` cut. Add
`--apply-y-max` only after that selection is intentionally frozen.

## Run and finalize

Run one prepared invocation:

```bash
python3 radiative_mode4.py run \
  mode4_rgk_pilot_s04468/manifest.json \
  --flat-index 4468 \
  --replica-index 0 \
  --executable build/aao_rad
```

After every prepared replica completes:

```bash
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
- core and unrestricted-tail fractions;
- explicit seed, event count, and `sigr_max`.

`run` writes products under a stem such as:

```text
runs/s04468/s04468__g0000
```

The `.norm` and run JSON record the stratum cross section, event weight,
proposal count, core/tail trial and event counts, maximum multiplicity, and
proposal efficiency. `aao_rad.mode4.csv` records accepted-event diagnostics
and is independently checked against the manifest proposal density and
final-coordinate bounds.

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

Compare `padding_0p03`, `padding_0p035`, and `padding_0p04`, initially at a
90/10 core/tail mixture. Inspect:

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
